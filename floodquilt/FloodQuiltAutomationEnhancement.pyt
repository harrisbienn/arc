# -*- coding: utf-8 -*-

"""
ArcGIS Pro Python Toolbox (.pyt) for Flood Quilt Automation
Filename: FloodQuiltAutomation.pyt
Author: Harris Bienn, Halff Associates, hbienn@halff.com
Date: 2025-09-26
Version: 0.1
--------------------------
Description:
This toolbox contains tools to automate the processing of flood hazard data,
including splitting by flood frequency, merging datasets, and creating topologies.
Each tool is designed to be modular, allowing for flexible workflows.
---------------------------
Requirements:
- ArcGIS Pro with arcpy
- Input data must be in a compatible format (feature classes within a geodatabase)
- Properly configured environment settings for geoprocessing
---------------------------
Notes:
- Ensure that the input feature classes have the expected schema, particularly the flood frequency field.
- The tools assume that the user has write access to the specified geodatabase.
- Error handling is minimal; users should validate inputs before running the tools.
- For best performance, run tools in an environment with sufficient resources.
- For inline variable substitution, parameters passed as a String are evaluated using locals(), globals() and isinstance(). To override, substitute values directly.
- To allow overwriting outputs change overwriteOutput option to True.
---------------------------
Warning:
*** Use arcpy.Parameter() with keyword args only (no positional args). Pro's validator is picky. ***
*** Use "Workspace" instead of "DEWorkspace" for GDB parameters in .pyt; some Pro installs reject DEWorkspace in Python toolboxes. ***
*** Use "Table View" instead of "DETable" for table parameters in .pyt; Table View is the most compatible datatype for user-picked tables in a .pyt. ***
*** Use arcpy.env to set environments, not param.env.
---------------------------
Changelog:
2025-09-26:     -Initial version, based on FloodQuiltSplitModule.py, FloodQuiltMergeModule.py, and FloodQuiltCombineModule.py
2025-09-29:     -Added Pipeline tool to run all steps in sequence. 
                -Enabled topology validation as an option.
2025-09-30:     -Refactored shared logic into helper functions. 
                -Improved parameter validation and messaging.
2025-10-01:     -Pairwise tools toggle (Dissolve/Erase)
                -Force 0-buffer pre-clean toggle
                -Keep intermediates toggle (handy for debugging)
                -Shared merge_floods() helper (so Tool 2 and Tool 4 share logic, like Split/Combine)
                -Central scratch FGDB reuse
                -Spatial index on heavy outputs
                -Timing decorator + run log so you get a mini perf report in GP messages
                -Wire-through of all toggles in Tool 2, 3, 4 (and Pipeline)
---------------------------
To Do:
- Add more robust error handling and input validation.
- Specifically topology error handing in the Merge submodule
- Allow user-defined field names for watershed and HUC8 in lookup table.
- Streamline temporary data handling to minimize I/O.
- 
"""

import arcpy
import os
import re
import time
import functools
from collections import defaultdict

# ---------------------------
# Global/Shared helpers
# ---------------------------

def _msg(msg, warn=False): 
    """Helper to print a message or warning."""
    (arcpy.AddWarning if warn else arcpy.AddMessage)(msg)

def choose_spatial_ref_from_inputs(paths, fallback_sr=None):
    """Choose spatial reference from first valid input path; else use WGS 1984 fallback."""
    if fallback_sr:
        return fallback_sr
    for p in paths:
        try:
            sr = arcpy.Describe(p).spatialReference
            if sr and sr.factoryCode not in (0, None):
                return sr
        except Exception:
            continue
    _msg("Could not detect spatial reference; defaulting to WGS 1984 (EPSG:4326).", warn=True)
    return arcpy.SpatialReference(4326)

def find_child_fc(workspace, candidates):
    """Find first existing FC in workspace from list of candidate names."""
    for name in candidates:
        candidate = os.path.join(workspace, name)
        if arcpy.Exists(candidate):
            return candidate
    return None

def apply_envs(parallel="100%", overwrite=True, workspace=None, scratch=None, ocs=None):
    """Apply common environment settings."""
    if overwrite is not None:
        arcpy.env.overwriteOutput = bool(overwrite)
    if parallel:
        # Accept "100%" or int
        arcpy.env.parallelProcessingFactor = parallel
    if workspace:
        arcpy.env.workspace = workspace
    if scratch:
        arcpy.env.scratchWorkspace = scratch
    if ocs:
        arcpy.env.outputCoordinateSystem = ocs

# ---------- Timing & logging ----------
def timed(label):
    def deco(fn):
        @functools.wraps(fn)
        def w(*a, **k):
            t0 = time.perf_counter()
            arcpy.AddMessage(f"[{label}] started")
            try:
                return fn(*a, **k)
            finally:
                dt = time.perf_counter() - t0
                arcpy.AddMessage(f"[{label}] completed in {dt:0.2f}s")
        return w
    return deco

class RunLog:
    def __init__(self): self.lines=[]
    def add(self, k, v): self.lines.append(f"{k}: {v}")
    def dump(self, header="Run summary"):
        arcpy.AddMessage(header)
        for ln in self.lines: arcpy.AddMessage("  " + ln)

# ---------- Scratch GDB manager (one per run) ----------
def get_scratch_gdb():
    sgdb = arcpy.env.scratchGDB
    if sgdb and arcpy.Exists(sgdb):
        return sgdb
    base_folder = (arcpy.env.scratchFolder or arcpy.env.workspace or arcpy.env.scratchWorkspace or os.getcwd())
    name = arcpy.CreateUniqueName("fgdb_", base_folder)
    if not name.lower().endswith(".gdb"):
        name = f"{name}.gdb"
    arcpy.management.CreateFileGDB(os.path.dirname(name), os.path.basename(name))
    return name

# ---------- Geometry helpers ----------
def add_spatial_index(fc):
    try:
        arcpy.management.AddSpatialIndex(fc)
    except Exception:
        pass

def _same_sr(fc, sr):
    try:
        dsr = arcpy.Describe(fc).spatialReference
        return (dsr and sr and dsr.factoryCode == sr.factoryCode and dsr.name == sr.name)
    except Exception:
        return False

def read_fc_paths_from_table(table, path_field): # returns list of strings
        paths = []
        with arcpy.da.SearchCursor(table, [path_field]) as cur:
            for (p,) in cur:
                if p:
                    paths.append(str(p))
        return paths

def group_by_suffix(paths): # returns dict: suffix -> list of paths
        SUFFIX_RE = re.compile(r'(?:_|)(500yr|100yr|10yr)$', re.IGNORECASE) # Matches _500yr, 500yr, _100yr, 100yr, _10yr, 10yr at end of string
        buckets = defaultdict(list)
        for p in paths:
            base = os.path.splitext(os.path.basename(p))[0]
            m = SUFFIX_RE.search(base)
            if m:
                key = m.group(1).lower()
                buckets[key].append(p)
            else:
                _msg(f"Skipped (no suffix match): {p}", warn=True)
        return buckets

def pick_spatial_ref(paths, fallback_sr=None): # returns arcpy.SpatialReference
        return choose_spatial_ref_from_inputs(paths, fallback_sr)

@timed("normalize_inputs")
def normalize_inputs(fc_list, target_sr, scratch_gdb, delete_null=True):
    out = []
    for i, fc in enumerate(fc_list, 1):
        if not arcpy.Exists(fc):
            continue
        tmp = arcpy.CreateUniqueName(f"norm_{i}", scratch_gdb)
        if _same_sr(fc, target_sr):
            arcpy.management.CopyFeatures(fc, tmp)
        else:
            arcpy.management.Project(fc, tmp, target_sr)
        try:
            arcpy.management.RepairGeometry(tmp, validation_method="OGC",
                                            delete_null="DELETE_NULL" if delete_null else "KEEP_NULL")
        except TypeError:
            arcpy.management.RepairGeometry(tmp, validation_method="OGC")
        out.append(tmp)
    return out

def preclean_polygon_fc(in_fc, scratch_gdb, force_buffer0=False):
    base = os.path.basename(in_fc)
    clean = arcpy.CreateUniqueName(f"clean_{base}", scratch_gdb)
    arcpy.management.CopyFeatures(in_fc, clean)
    try:
        arcpy.management.RepairGeometry(clean, validation_method="OGC", delete_null="DELETE_NULL")
    except TypeError:
        arcpy.management.RepairGeometry(clean, validation_method="OGC")
    gtype = arcpy.Describe(clean).shapeType
    if force_buffer0 and gtype and gtype.lower() == "polygon":
        buf = arcpy.CreateUniqueName(f"buf0_{base}", scratch_gdb)
        arcpy.analysis.Buffer(clean, buf, "0 Feet", dissolve_option="NONE")
        return buf
    return clean

# ---------- Pairwise-capable overlays ----------
def pairwise_available():
    try:
        # Simple capability probe
        _ = arcpy.analysis.PairwiseDissolve
        _ = arcpy.analysis.PairwiseErase
        return True
    except Exception:
        return False

@timed("fast_dissolve")
def fast_dissolve(in_fc, out_fc, field, use_pairwise=True, is_polyline=False):
    used_pairwise = False
    if use_pairwise and pairwise_available():
        try:
            arcpy.analysis.PairwiseDissolve(in_fc, out_fc, field)
            used_pairwise = True
        except Exception:
            used_pairwise = False
    if not used_pairwise:
        kwargs = dict(in_features=in_fc, out_feature_class=out_fc,
                    dissolve_field=field, statistics_fields=[], multi_part="MULTI_PART")
        if is_polyline:
            kwargs["unsplit_lines"] = "DISSOLVE_LINES"
        arcpy.management.Dissolve(**kwargs)

@timed("safe_dissolve")
def safe_dissolve(in_fc, out_fc, field, use_pairwise=True, force_buffer0=False):
    gtype = (arcpy.Describe(in_fc).shapeType or "").lower()
    is_polyline = (gtype == "polyline")
    try:
        fast_dissolve(in_fc, out_fc, field, use_pairwise=use_pairwise, is_polyline=is_polyline)
        return out_fc
    except arcpy.ExecuteError:
        msg = arcpy.GetMessages()
        if "160196" in msg or "Invalid Topology" in msg:
            arcpy.AddWarning("Dissolve failed (160196). Pre-cleaning + retry…")
            clean = preclean_polygon_fc(in_fc, get_scratch_gdb(), force_buffer0=force_buffer0)
            fast_dissolve(clean, out_fc, field, use_pairwise=use_pairwise, is_polyline=is_polyline)
            return out_fc
        raise

@timed("fast_erase")
def fast_erase(in_fc, erase_fc, out_fc, use_pairwise=True):
    used_pairwise = False
    if use_pairwise and pairwise_available():
        try:
            arcpy.analysis.PairwiseErase(in_fc, erase_fc, out_fc)
            used_pairwise = True
        except Exception:
            used_pairwise = False
    if not used_pairwise:
        arcpy.analysis.Erase(in_fc, erase_fc, out_fc)

# ---------- Basename ----------

def sanitize_alnum(s: str) -> str:
    """Strict sanitize: keep letters and numbers only (no spaces/punct)."""
    return re.sub(r'[^A-Za-z0-9]', '', s or '')

def compute_basename(watershed: str, huc8: str) -> str:
    """Basename = <Watershed>_<HUC8> (both sanitized)."""
    return f"{sanitize_alnum(watershed)}_{sanitize_alnum(huc8)}"

''' # Future enhancement: allow user to specify field names
def build_ws_huc_lookup(table_path: str, ws_field: str, huc_field: str): 
    """Return a dictionary mapping watershed to a sorted list of HUC8s."""
    ws_to_huc = defaultdict(set)

    with arcpy.da.SearchCursor(table_path, [ws_field, huc_field]) as cur:
        for w, h in cur:
            if w is not None and h is not None:
                ws_to_huc[str(w)].add(str(h))

    return {k: sorted(list(v)) for k, v in ws_to_huc.items()}
'''

def build_ws_huc_lookup(table_path: str, ws_field_guess="Watershed", huc_field_guess="Name"): # Most common failure point is field names not matching guesses
    """Return dict: watershed -> sorted list of HUC8s."""
    ws_to_huc = defaultdict(set)
    fields = [f.name for f in arcpy.ListFields(table_path)]

    def findcol(want):
        """Find column name in fields, case-insensitive match; else return want."""
        for f in fields:
            if f.lower() == want.lower():
                return f
        return want

    ws_col  = findcol(ws_field_guess)
    huc_col = findcol(huc_field_guess)

    with arcpy.da.SearchCursor(table_path, [ws_col, huc_col]) as cur:
        for w, h in cur:
            if w is not None and h is not None:
                ws_to_huc[str(w)].add(str(h))
    return {k: sorted(list(v)) for k, v in ws_to_huc.items()}

def _updateWatershedAndHUC8Filters(p_lut, p_watershed, p_huc8):
    """ Updates the filter lists for Watershed and HUC8 based on the Lookup Table """
    lookup = build_ws_huc_lookup(p_lut.valueAsText)
    p_watershed.filter.list = sorted(lookup.keys())
    if p_watershed.value and p_watershed.value in lookup:
        p_huc8.filter.list = lookup[p_watershed.value]
    else:
        p_huc8.filter.list = []

def _updateDerivedPaths(p_gdb, p_watershed, p_huc8, p_basename, p_outfds, p_flood_500, p_flood_100, p_flood_10):
    """ Computes and assigns derived paths for output datasets based on watershed and HUC8 selections """
    basename_value = compute_basename(p_watershed.value, p_huc8.value)
    p_basename.value = basename_value
    base_path = os.path.join(p_gdb.valueAsText, f"Output_{basename_value}")
    p_outfds.value = base_path
    p_flood_500.value = os.path.join(base_path, f"Existing_Future_{basename_value}_500yr")
    p_flood_100.value = os.path.join(base_path, f"Existing_Future_{basename_value}_100yr")
    p_flood_10.value = os.path.join(base_path, f"Existing_Future_{basename_value}_10yr")

# -----------Module helpers-----------
@timed("split_floods")
def split_floods(gdb: str, 
                 basename: str, 
                 spatial_ref, 
                 existing_fc: str, 
                 future_fc: str, 
                 split_field: str = "FLOOD_FREQ"):
    """
    Core Split step shared by Tool 1 (SplitFloodFreqs) and Tool 4 (Pipeline).

    Creates feature datasets:
      Existing_<basename>, Future_<basename>
    Runs SplitByAttributes on <split_field>, converts likely parts to single-part, and deletes
    the raw split outputs. Returns a dict with paths to datasets and single-part outputs.

    Returns:
      {
        "fds_exist": <path>,
        "fds_future": <path>,
        "exist_10": <path or "">,
        "exist_100": <path or "">,
        "exist_500": <path or "">,
        "future_10": <path or "">,
        "future_100": <path or "">,
        "future_500": <path or "">
      }
    """
    # Ensure datasets
    fds_exist_name = f"Existing_{basename}"
    fds_future_name = f"Future_{basename}"
    fds_exist = os.path.join(gdb, fds_exist_name)
    fds_future = os.path.join(gdb, fds_future_name)

    if not arcpy.Exists(fds_exist):
        _msg(f"Creating feature dataset: {fds_exist}")
        arcpy.management.CreateFeatureDataset(gdb, fds_exist_name, spatial_ref)
    if not arcpy.Exists(fds_future):
        _msg(f"Creating feature dataset: {fds_future}")
        arcpy.management.CreateFeatureDataset(gdb, fds_future_name, spatial_ref)

    # ---- Process the EXISTING feature class ----
    # Split by attributes
    _msg(f"Splitting existing feature class by {split_field} and converting multipart to singlepart…")
    arcpy.analysis.SplitByAttributes(existing_fc, fds_exist, [split_field])

    # Helper to find the split children regardless of how SplitByAttributes named them
    def _pick(ws, candidates):
        for c in candidates:
            p = os.path.join(ws, c)
            if arcpy.Exists(p):
                return p
        return None

    # Support both T10/T1/T0_2 and 10yr/100yr/500yr (and a few variants)
    E10 = _pick(fds_exist, ["T10", "10yr", "10", "YR10", "yr10"])
    E1 = _pick(fds_exist, ["T1", "100yr", "100", "YR100", "yr100"])
    E5 = _pick(fds_exist, ["T0_2", "500yr", "500", "YR500", "yr500"])

    # Canonical single-part outputs for existing
    exist_10 = os.path.join(fds_exist, f"Existing_{basename}_10yr") if E10 else ""
    exist_100 = os.path.join(fds_exist, f"Existing_{basename}_100yr") if E1 else ""
    exist_500 = os.path.join(fds_exist, f"Existing_{basename}_500yr") if E5 else ""

    # Multipart → Singlepart for existing
    for src, dst in ((E10, exist_10), (E1, exist_100), (E5, exist_500)):
        if src and dst:
            arcpy.management.MultipartToSinglepart(src, dst)

    # Delete raw split outputs if converted for existing
    to_delete_existing = [p for p in (E10, E1, E5) if p and arcpy.Exists(p)]
    if to_delete_existing:
        arcpy.management.Delete(to_delete_existing)
    _msg("Existing feature class processing complete.")

    # ---- Process the FUTURE feature class ----
    # Split by attributes
    _msg(f"Splitting future feature class by {split_field} and converting multipart to singlepart…")
    arcpy.analysis.SplitByAttributes(future_fc, fds_future, [split_field])

    F10 = _pick(fds_future, ["T10", "10yr", "10", "YR10", "yr10"])
    F1 = _pick(fds_future, ["T1", "100yr", "100", "YR100", "yr100"])
    F5 = _pick(fds_future, ["T0_2", "500yr", "500", "YR500", "yr500"])

    # Canonical single-part outputs for future
    future_10 = os.path.join(fds_future, f"Future_{basename}_10yr") if F10 else ""
    future_100 = os.path.join(fds_future, f"Future_{basename}_100yr") if F1 else ""
    future_500 = os.path.join(fds_future, f"Future_{basename}_500yr") if F5 else ""

    # Multipart → Singlepart for future
    for src, dst in ((F10, future_10), (F1, future_100), (F5, future_500)):
        if src and dst:
            arcpy.management.MultipartToSinglepart(src, dst)

    # Delete raw split outputs if converted for future
    to_delete_future = [p for p in (F10, F1, F5) if p and arcpy.Exists(p)]
    if to_delete_future:
        arcpy.management.Delete(to_delete_future)
    _msg("Future feature class processing complete.\nSplit step complete.")

    _msg(f"Outputs:\n  Geodatabase: {gdb}\n  Existing 500yr: {exist_500}\n  Existing 100yr: {exist_100}\n  Existing 10yr: {exist_10}\n  Future 500yr: {future_500}\n  Future 100yr: {future_100}\n  Future 10yr: {future_10}")
    
    return {
        "fds_exist": fds_exist,
        "fds_future": fds_future,
        "exist_10": exist_10 if (exist_10 and arcpy.Exists(exist_10)) else "",
        "exist_100": exist_100 if (exist_100 and arcpy.Exists(exist_100)) else "",
        "exist_500": exist_500 if (exist_500 and arcpy.Exists(exist_500)) else "",
        "future_10": future_10 if (future_10 and arcpy.Exists(future_10)) else "",
        "future_100": future_100 if (future_100 and arcpy.Exists(future_100)) else "",
        "future_500": future_500 if (future_500 and arcpy.Exists(future_500)) else "",
    }

@timed("merge_floods")
def merge_floods(
    project_gdb: str,
    basename: str,
    buckets: dict,               # {"500yr":[paths], "100yr":[paths], "10yr":[paths]}
    target_sr,
    use_pairwise: bool = False,
    force_buffer0: bool = False,
    keep_intermediates: bool = False
):
    scratch = get_scratch_gdb()

    # Create Output_<Basename>
    fds_name = f"Output_{basename}"
    out_fds = os.path.join(project_gdb, fds_name)
    if not arcpy.Exists(out_fds):
        _msg(f"Creating feature dataset: {out_fds}")
        arcpy.management.CreateFeatureDataset(project_gdb, fds_name, target_sr)

    '''# Temp merged paths
    merged_500 = os.path.join(out_fds, "Merged_500yr_tmp")
    merged_100 = os.path.join(out_fds, "Merged_100yr_tmp")
    merged_010 = os.path.join(out_fds, "Merged_10yr_tmp")

    # Final outputs (Keep classic names to match Combine)
    out_500 = os.path.join(out_fds, "Flood_500yr")
    out_100 = os.path.join(out_fds, "Flood_100yr")
    out_010 = os.path.join(out_fds, "Flood_10yr")'''

    # Merge temps
    merged_500 = os.path.join(out_fds, f"Merged_{basename}_500yr")
    merged_100 = os.path.join(out_fds, f"Merged_{basename}_100yr")
    merged_010 = os.path.join(out_fds, f"Merged_{basename}_10yr")

    # Final
    out_500 = os.path.join(out_fds, f"Existing_Future_{basename}_500yr")
    out_100 = os.path.join(out_fds, f"Existing_Future_{basename}_100yr")
    out_010 = os.path.join(out_fds, f"Existing_Future_{basename}_10yr")

    def _merge_one(key, merged_path, out_path):
        if not buckets.get(key):
            return None, None
        norm = normalize_inputs(buckets[key], target_sr, scratch_gdb=scratch, delete_null=True)
        if not norm:
            return None, None
        _msg(f"Merging {key} ({len(norm)} inputs)…")
        arcpy.management.Merge(norm, merged_path)
        add_spatial_index(merged_path)
        
        # Safe dissolve (pairwise if available)
        _msg(f"Dissolving merged {key} on FLOOD_FREQ (pairwise={use_pairwise}, force_buffer0={force_buffer0})…")
        safe_dissolve(merged_path, out_path, "FLOOD_FREQ",
                      use_pairwise=use_pairwise, force_buffer0=force_buffer0)
        add_spatial_index(out_path)
        return merged_path, out_path

    m500, o500 = _merge_one("500yr", merged_500, out_500)
    m100, o100 = _merge_one("100yr", merged_100, out_100)
    m010, o010 = _merge_one("10yr",  merged_010, out_010)

    if not keep_intermediates:
        for tmp in (m500, m100, m010):
            if tmp and arcpy.Exists(tmp):
                try: arcpy.management.Delete(tmp)
                except Exception: pass

    _msg(f"""Merge step complete.\n Outputs:
         \n     Merged 500yr: {m500 if (m500 and arcpy.Exists(m500)) else 'Intermediate output not preserved'}
         \n     Merged 100yr: {m100 if (m100 and arcpy.Exists(m100)) else 'Intermediate output not preserved'}
         \n     Merged 10yr: {m010 if (m010 and arcpy.Exists(m010)) else 'Intermediate output not preserved'}
         \n     Dissolved 500yr: {o500 if (o500 and arcpy.Exists(o500)) else 'N/A'}
         \n     Dissolved 100yr: {o100 if (o100 and arcpy.Exists(o100)) else 'N/A'}
         \n     Dissolved 10yr: {o010 if (o010 and arcpy.Exists(o010)) else 'N/A'}""")
    
    return {"out_fds": out_fds,
            "out_500": o500 if (o500 and arcpy.Exists(o500)) else "",
            "out_100": o100 if (o100 and arcpy.Exists(o100)) else "",
            "out_010": o010 if (o010 and arcpy.Exists(o010)) else ""}

@timed("combine_floods")
def combine_floods(
    gdb: str,
    basename: str,
    out_fds: str,   # Output_<Basename>
    fc_500: str,
    fc_100: str,
    fc_10: str,
    out_cs=None,
    validate_topology: bool = True,
    use_pairwise: bool = True,
    force_buffer0: bool = False,
    keep_intermediates: bool = False
    ):
    scratch = get_scratch_gdb()
    
    # Final dataset (SR from inputs or env/param)
    final_fds_name = f"Final_{basename}"
    final_fds = os.path.join(gdb, final_fds_name)
    sr_final = None
    for fc in (fc_10, fc_100, fc_500):
        if fc and arcpy.Exists(fc):
            dsr = arcpy.Describe(fc).spatialReference
            if dsr and dsr.factoryCode:
                sr_final = dsr; break
    if sr_final is None:
        _msg("Could not determine spatial reference from inputs; using environment setting.", warn=True)
        sr_final = out_cs or arcpy.env.outputCoordinateSystem
    if not arcpy.Exists(final_fds):
        _msg(f"Creating feature dataset: {final_fds} with spatial reference: {sr_final.name if sr_final else 'Unknown'}")
        arcpy.management.CreateFeatureDataset(gdb, final_fds_name, sr_final)

    # Repair inputs & (optional) 0-buffer
    _msg(f"Pre-cleaning input feature classes…")
    rep_500 = preclean_polygon_fc(fc_500, scratch, force_buffer0=force_buffer0) if (fc_500 and arcpy.Exists(fc_500)) else None
    rep_100 = preclean_polygon_fc(fc_100, scratch, force_buffer0=force_buffer0) if (fc_100 and arcpy.Exists(fc_100)) else None
    rep_010 = preclean_polygon_fc(fc_10,  scratch, force_buffer0=force_buffer0) if (fc_10  and arcpy.Exists(fc_10))  else None

    # Erase overlaps (store in Output_<Basename>)
    erased_500 = os.path.join(out_fds, f"{basename}_500yr_Erased")
    erased_100 = os.path.join(out_fds, f"{basename}_100yr_Erased")
    erased_10 = os.path.join(out_fds, f"{basename}_10yr_Erased")
    
    _msg("Erasing 500yr - 100yr...")
    if rep_500 and rep_100: fast_erase(rep_500, rep_100, erased_500, use_pairwise=use_pairwise)
    else: erased_500 = ""
    _msg("Erasing 100yr - 10yr...")
    if rep_100 and rep_010: fast_erase(rep_100, rep_010, erased_100, use_pairwise=use_pairwise)
    else: erased_100 = ""

    # Merge the erased 500yr and 100yr
    ef_500_100 = os.path.join(out_fds, f"{basename}_500yr_100yr_Erased")
    arcpy.managment.Merge([erased_500, erased_100], ef_500_100)

    # Erase the merged 500yr and 100yr from the 10yr 
    if ef_500_100 and rep_010: fast_erase(rep_010, ef_500_100, erased_10, use_pairwise=use_pairwise)

    # Merge final thematic stack
    ef_merged = os.path.join(final_fds, f"Existing_Future_{basename}_500yr_100yr_10yr")
    inputs = [fc for fc in (rep_010, erased_500 if erased_500 else None, erased_100 if erased_100 else None) if fc and arcpy.Exists(fc)]
    if not inputs:
        raise RuntimeError("Combine step: nothing to merge. Ensure Merge outputs exist.")
    _msg(f"Merging final thematic stack ({len(inputs)} inputs)…")
    arcpy.management.Merge(inputs, ef_merged)
    add_spatial_index(ef_merged)

    # Topology build (rule always added)
    _msg("Creating topology and adding 'Must Not Overlap (Area)' rule…")
    topo = arcpy.management.CreateTopology(final_fds, f"Final_{basename}_Topology").getOutput(0)
    arcpy.management.AddFeatureClassToTopology(topo, ef_merged, 1, 1)
    arcpy.management.AddRuleToTopology(topo, "Must Not Overlap (Area)", ef_merged)
    if validate_topology:
        arcpy.AddMessage("Validating topology (Full Extent)…")
        arcpy.management.ValidateTopology(topo, "Full_Extent")

    # Cleanup
    if not keep_intermediates:
        _msg("Cleaning up intermediate datasets…")
        for fc in (rep_500, rep_100, rep_010, erased_500 if erased_500 else None, erased_100 if erased_100 else None):
            if fc and arcpy.Exists(fc): 
                try: arcpy.management.Delete(fc)
                except Exception: pass
    
    _msg(f"""Combine step complete.\n Outputs:
         \n     Erased 500yr-100yr: {erased_500 if (erased_500 and arcpy.Exists(erased_500)) else 'Intermediate output not preserved'}
         \n     Erased 100yr-10yr: {erased_100 if (erased_100 and arcpy.Exists(erased_100)) else 'Intermediate output not preserved'}
         \n     Final merged 500yr/100yr/10yr: {ef_merged if (ef_merged and arcpy.Exists(ef_merged)) else 'N/A'}
         \n     Topology: {topo if (topo and arcpy.Exists(topo)) else 'N/A'}""")  

    return {"final_fds": final_fds, "topology": topo, "ef_merged": ef_merged,
            "erased_500": erased_500, "erased_100": erased_100}


# ---------------------------
# Toolbox container (tool order matters)
# ---------------------------
class Toolbox(object):
    def __init__(self):
        self.label = "Flood Quilt Automation"
        self.alias = "Flood Quilt Automation"
        """1) Split, 2) Merge, 3) Combine, 4) Pipeline"""
        self.tools = [SplitFloodFreqs, MergeFloodFreqs, CombineFloodFreqs, RunFloodQuiltPipeline]

# ---------------------------
# Tool 1: SplitFloodFreqs based on FloodQuiltSplitModule.py
# ---------------------------
class SplitFloodFreqs(object):
    """
    Create Existing_<Basename> and Future_<Basename> FDS in Project GDB,
    SplitByAttributes on FLOOD_FREQ (or chosen field) 
    Single-part outputs:
      Existing_<Basename>_{10yr,100yr,500yr}
      Future_<Basename>_{10yr,100yr,500yr}
    Derived output: Basename
    """
    def __init__(self):
        self.label = "1) Split Flood Frequencies"
        self.description = "Split Existing/Future flood layers by frequency into Existing/Future datasets with 10yr/100yr/500yr outputs."
        self.canRunInBackground = False

    def getParameterInfo(self):
        # Selector
        p_lut = arcpy.Parameter(
            displayName="Lookup Table (two columns: Watershed, HUC8)",
            name="lookup_table",
            datatype="Table View",          # ← was DETable, Table View is the most compatible datatype for user-picked tables in a .pyt.
            parameterType="Required",
            direction="Input"
        )
        p_watershed  = arcpy.Parameter(
            displayName="Watershed",
            name="watershed",
            datatype="GPString",
            parameterType="Required",
            direction="Input"
        ); p_watershed.filter.type = "ValueList"
        p_huc8 = arcpy.Parameter(
            displayName="HUC8",
            name="huc8",
            datatype="GPString",
            parameterType="Required",
            direction="Input"
        ); p_huc8.filter.type = "ValueList"

        # Processing
        p_gdb = arcpy.Parameter(
            displayName="Project Geodatabase",
            name="project_gdb",
            datatype="Workspace",           # ← was DEWorkspace, Workspace is the safest for GDBs in a .pyt.
            parameterType="Required",
            direction="Input"
        )
        p_sr  = arcpy.Parameter(
            displayName="Target Spatial Reference",
            name="target_sr",
            datatype="GPSpatialReference",
            parameterType="Required",
            direction="Input"
        )
        p_in_exist  = arcpy.Parameter(
            displayName="Existing Flood Feature Class",
            name="input_feature_existing",
            datatype="DEFeatureClass",
            parameterType="Required",
            direction="Input"
        )
        p_in_future = arcpy.Parameter(
            displayName="Future Flood Feature Class",
            name="input_feature_future",
            datatype="DEFeatureClass",
            parameterType="Required",
            direction="Input"
        )
        p_splitfield= arcpy.Parameter(
            displayName="Split Field (frequency)",
            name="split_field",
            datatype="GPString",
            parameterType="Optional",
            direction="Input"
        ); p_splitfield.value = "FLOOD_FREQ"

        # Optional envs
        p_parallel = arcpy.Parameter(
            displayName="Parallel Processing Factor",
            name="parallel_factor",
            datatype="GPString",
            parameterType="Optional",
            direction="Input",
            category="Environment Settings"
        ); p_parallel.value = "100%"
        p_overwrite= arcpy.Parameter(       # Check functionality, cascade through to applicable tools
            displayName="Overwrite Outputs",
            name="overwrite_outputs",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input",
            category="Environment Settings"
        ); p_overwrite.value = True

        # Derived
        p_basename = arcpy.Parameter(
            displayName="Basename (derived)",
            name="basename",
            datatype="GPString",
            parameterType="Derived",
            direction="Output"
        )
        p_fds_exist = arcpy.Parameter(
            displayName="Existing Feature Dataset",
            name="existing_feature_dataset",
            datatype="DEFeatureDataset",
            parameterType="Derived",
            direction="Output"
        )
        p_exist_10  = arcpy.Parameter(
            displayName="Existing 10yr",
            name="existing_10yr",
            datatype="DEFeatureClass",
            parameterType="Derived",
            direction="Output"
        )
        p_exist_100 = arcpy.Parameter(
            displayName="Existing 100yr",
            name="existing_100yr",
            datatype="DEFeatureClass",
            parameterType="Derived",
            direction="Output"
        )
        p_exist_500 = arcpy.Parameter(
            displayName="Existing 500yr",
            name="existing_500yr",
            datatype="DEFeatureClass",
            parameterType="Derived",
            direction="Output"
        )
        p_fds_future= arcpy.Parameter(
            displayName="Future Feature Dataset",
            name="future_feature_dataset",
            datatype="DEFeatureDataset",
            parameterType="Derived",
            direction="Output"
        )
        p_future_10 = arcpy.Parameter(
            displayName="Future 10yr",
            name="future_10yr",
            datatype="DEFeatureClass",
            parameterType="Derived",
            direction="Output"
        )
        p_future_100= arcpy.Parameter(
            displayName="Future 100yr",
            name="future_100yr",
            datatype="DEFeatureClass",
            parameterType="Derived",
            direction="Output"
        )
        p_future_500= arcpy.Parameter(
            displayName="Future 500yr",
            name="future_500yr",
            datatype="DEFeatureClass",
            parameterType="Derived",
            direction="Output"
        )

        return [p_lut, p_watershed, p_huc8, p_gdb, p_sr, p_in_exist, p_in_future, p_splitfield,
            p_parallel, p_overwrite, 
            p_basename, p_fds_exist, p_exist_10, p_exist_100, p_exist_500,
            p_fds_future, p_future_10, p_future_100, p_future_500]

    def updateParameters(self, params):
        p_lut = params[0]
        p_watershed = params[1]
        p_huc8 = params[2]

        # Cascade updates based on lookup table
        if p_lut.altered and p_lut.value:
            _updateWatershedAndHUC8Filters(p_lut, p_watershed, p_huc8)
        return

    def execute(self, params, messages):
        (p_lut, p_watershed, p_huc8, p_gdb, p_sr, p_in_exist, p_in_future, p_splitfield,
         p_parallel, p_overwrite, 
         p_basename, p_fds_exist, p_exist_10, p_exist_100, p_exist_500,
         p_fds_future, p_future_10, p_future_100, p_future_500) = params

        apply_envs(p_parallel.valueAsText or "100%",
                   p_overwrite.value,  
                   )

        gdb         = p_gdb.valueAsText
        spatial_ref = p_sr.value
        in_exist    = p_in_exist.valueAsText
        in_future   = p_in_future.valueAsText
        split_field = (p_splitfield.valueAsText or "FLOOD_FREQ")
        
        # Basename: prefer explicit; else derive from WS/HUC8
        basename = p_basename.valueAsText or (compute_basename(p_watershed.valueAsText, p_huc8.valueAsText) if (p_watershed.value and p_huc8.value) else None)
        if not basename:
            raise RuntimeError("Basename is required (either supply it or choose Watershed + HUC8).")
        p_basename.value = basename

        # Shared core
        r = split_floods(gdb, 
                        basename, 
                        spatial_ref, 
                        in_exist, 
                        in_future, 
                        split_field
                        )
        
        # Push derived outputs
        p_fds_exist.value  = r["fds_exist"]
        p_fds_future.value = r["fds_future"]
        p_exist_10.value   = r["exist_10"]
        p_exist_100.value  = r["exist_100"]
        p_exist_500.value  = r["exist_500"]
        p_future_10.value  = r["future_10"]
        p_future_100.value = r["future_100"]
        p_future_500.value = r["future_500"]
        
        
        # Paths to feature classes that were modified
        fc_paths = {
            "Existing_10yr": r["exist_10"],
            "Existing_100yr": r["exist_100"],
            "Existing_500yr": r["exist_500"],
            "Future_10yr": r["future_10"],
            "Future_100yr": r["future_100"],
            "Future_500yr": r["future_500"],
        }

        # Path for output table
        output_table = os.path.join(gdb, "FeatureClassPaths")
        
        # Create table if not exists
        arcpy.management.CreateTable(gdb, "FeatureClassPaths")
        arcpy.management.AddField(output_table, "Name", "TEXT")
        arcpy.management.AddField(output_table, "FCPath", "TEXT")
        
        # Insert paths into table
        with arcpy.da.InsertCursor(output_table, ["Name", "FCPath"]) as cursor:
            for name, path in fc_paths.items():
                cursor.insertRow([name, path])

        _msg(f"Feature class paths have been recorded in {output_table}.")

# ---------------------------
# Tool 2: MergeFloodFreqs based on FloodQuiltMergeModule.py
# ---------------------------
class MergeFloodFreqs(object):
    """
    Read table of FC paths, group by suffix (500yr, 100yr, 10yr), merge each group, 
    dissolve by FLOOD_FREQ into Output_<Basename> in Project GDB.
    Inputs: table with FC paths, Project GDB, Basename (from Split tool)
    Derived outputs: Output_500yr, Output_100yr, Output_10yr
    """
    SUFFIX_RE = re.compile(r'(?:_|)(500yr|100yr|10yr)$', re.IGNORECASE) # Matches _500yr, 500yr, _100yr, 100yr, _10yr, 10yr at end of string
    FREQ_KEYS = ["500yr", "100yr", "10yr"] # order matters for output naming
    DISSOLVE_FIELD = "FLOOD_FREQ" # Field to dissolve on; ensure this field exists in inputs

    def __init__(self):
        self.label = "2) Merge Flood Frequencies"
        self.description = "Merge 6 inputs into 3 by suffix and dissolve by FLOOD_FREQ into Output feature dataset."
        self.canRunInBackground = False

    def getParameterInfo(self):
        p_lut = arcpy.Parameter(
            displayName="Lookup Table (two columns: Watershed, HUC8)",
            name="lookup_table",
            datatype="Table View",          # ← was DETable, Table View is the most compatible datatype for user-picked tables in a .pyt.
            parameterType="Required",
            direction="Input"
        )
        p_watershed  = arcpy.Parameter(
            displayName="Watershed",
            name="watershed",
            datatype="GPString",
            parameterType="Required",
            direction="Input"
        ); p_watershed.filter.type = "ValueList"
        p_huc8 = arcpy.Parameter(
            displayName="HUC8",
            name="huc8",
            datatype="GPString",
            parameterType="Required",
            direction="Input"
        ); p_huc8.filter.type = "ValueList"
        p_table = arcpy.Parameter(
            displayName="Input Table (with FC paths)",
            name="input_table",
            datatype="Table View",          # ← was DETable, Table View is the most compatible datatype for user-picked tables in a .pyt.
            parameterType="Required",
            direction="Input"
        )
        p_field = arcpy.Parameter(
            displayName="Path Field Name",
            name="path_field",
            datatype="GPString",
            parameterType="Optional",
            direction="Input"
        ); p_field.value = "FCPath"
        p_gdb   = arcpy.Parameter(
            displayName="Project Geodatabase",
            name="project_gdb",
            datatype="Workspace",           # ← was DEWorkspace, Workspace is the safest for GDBs in a .pyt.
            parameterType="Required",
            direction="Input"
        )
        p_basename  = arcpy.Parameter(         # Check if being passed correctly
            displayName="Basename (from Split tool)",
            name="basename",
            datatype="GPString",
            parameterType="Derived",
            direction="Input"
        )
        p_sr    = arcpy.Parameter(
            displayName="Target Spatial Reference (optional)",
            name="target_sr",
            datatype="GPSpatialReference",
            parameterType="Optional",
            direction="Input"
        )

        # Optional envs
        p_parallel = arcpy.Parameter(
            displayName="Parallel Processing Factor",
            name="parallel_factor",
            datatype="GPString",
            parameterType="Optional",
            direction="Input",
            category="Environment Settings"
        ); p_parallel.value = "100%"
        p_overwrite= arcpy.Parameter(      # Check functionality, cascade through to applicable tools 
            displayName="Overwrite Outputs",
            name="overwrite_outputs",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input",
            category="Environment Settings"
        ); p_overwrite.value = True
        p_use_pairwise = arcpy.Parameter(
            displayName="Use pairwise tools where possible.\n(Faster option but may fail on complex geometries.)",
            name="use_pairwise",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input",
            category="Debug"
        ); p_use_pairwise.value = False
        p_force_buf0 = arcpy.Parameter(
            displayName="Force 0-buffer pre-clean before erase/merge.\n(Slower option but topologically safer.)",
            name="force_buffer0",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input",
            category="Debug"
        ); p_force_buf0.value = False
        p_keep_inter = arcpy.Parameter(
            displayName="Keep intermediate datasets",
            name="keep_intermediates",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input",
            category="Debug"
        ); p_keep_inter.value = False

        # Derived
        p_outfds = arcpy.Parameter(
            displayName="Output Feature Dataset",
            name="out_fds",
            datatype="DEFeatureDataset",
            parameterType="Derived",
            direction="Output"
        )
        p_out500 = arcpy.Parameter(
            displayName="Flood_500yr",
            name="out_500yr",
            datatype="DEFeatureClass",
            parameterType="Derived",
            direction="Output"
        )
        p_out100 = arcpy.Parameter(
            displayName="Flood_100yr",
            name="out_100yr",
            datatype="DEFeatureClass",
            parameterType="Derived",
            direction="Output"
        )
        p_out10  = arcpy.Parameter(
            displayName="Flood_10yr",
            name="out_10yr",
            datatype="DEFeatureClass",
            parameterType="Derived",
            direction="Output"
        )

        return [p_lut, p_watershed, p_huc8, p_table, p_field, p_gdb, p_basename, p_sr,
                p_parallel, p_overwrite, 
                p_use_pairwise, p_force_buf0, p_keep_inter,
                p_outfds, p_out500, p_out100, p_out10]

    def updateParameters(self, params):
        p_lut = params[0]
        p_watershed = params[1]
        p_huc8= params[2] 
        p_basename = params[6]

        # Cascade
        if p_lut.altered and p_lut.value:
            _updateWatershedAndHUC8Filters(p_lut, p_watershed, p_huc8)
        return

    def execute(self, params, messages):
        (p_lut, p_watershed, p_huc8, p_table, p_field, p_gdb, p_basename, p_sr,
        p_parallel, p_overwrite,
        p_use_pairwise, p_force_buf0, p_keep_inter,
        p_outfds, p_out500, p_out100, p_out10) = params

        apply_envs(p_parallel.valueAsText or "100%", p_overwrite.value)

        # Basename: prefer explicit; else derive from WS/HUC8
        basename = p_basename.valueAsText or (compute_basename(p_watershed.valueAsText, p_huc8.valueAsText) if (p_watershed.value and p_huc8.value) else None)
        if not basename:
            raise RuntimeError("Basename is required (either supply it or choose Watershed + HUC8).")
        p_basename.value = basename

        # Read + group paths
        paths = read_fc_paths_from_table(p_table.valueAsText, p_field.valueAsText or "FCPath")
        if len(paths) < 3:
            raise RuntimeError("Expected at least 3 rows with FC paths.")
        buckets = group_by_suffix(paths)
        for key in self.FREQ_KEYS:
            if not buckets.get(key): _msg(f"No inputs found for '{key}' group.", warn=True)

        sr = pick_spatial_ref(paths, p_sr.value)

        r = merge_floods(project_gdb=p_gdb.valueAsText,
                        basename=basename,
                        buckets=buckets,
                        target_sr=sr,
                        use_pairwise=bool(p_use_pairwise.value),
                        force_buffer0=bool(p_force_buf0.value),
                        keep_intermediates=bool(p_keep_inter.value))

        # Derived outputs
        p_outfds.value  = r["out_fds"]
        p_out500.value  = r["out_500"]
        p_out100.value  = r["out_100"]
        p_out10.value   = r["out_010"]

# ---------------------------
# Tool 3: CombineFloodFreqs ← based on FloodQuiltCombineModule.py
# ---------------------------
class CombineFloodFreqs(object):
    """
    Final combination step:
      - Create Final_<Basename> FDS
      - Repair geometry on Flood_500yr/100yr/10yr
      - Erase: 500yr-100yr, 100yr-10yr
      - Merge: (10yr, 500yr_erased, 100yr_erased) → Existing_Future_<Basename>_500yr_100yr_10yr
      - Create topology and add 'Must Not Overlap (Area)' rule
    Inputs: Project GDB, Basename (from Split tool), Output Feature Dataset (from Merge tool), Flood_500yr, Flood_100yr, Flood_10yr
    Derived outputs: Final Feature Dataset, Final Topology, Existing_Future_<Basename>_500yr_100yr_10yr
    """
    def __init__(self):
        self.label = "3) Combine Flood Frequencies"
        self.description = "Repairs, erases overlaps, merges, and builds topology for final features."
        self.canRunInBackground = False

    def getParameterInfo(self):
        p_lut = arcpy.Parameter(
            displayName="Lookup Table (two columns: Watershed, HUC8)",
            name="lookup_table",
            datatype="Table View",          # ← was DETable, Table View is the most compatible datatype for user-picked tables in a .pyt.
            parameterType="Required",
            direction="Input"
        )
        p_watershed  = arcpy.Parameter(
            displayName="Watershed",
            name="watershed",
            datatype="GPString",
            parameterType="Required",
            direction="Input"
        ); p_watershed.filter.type = "ValueList"
        p_huc8  = arcpy.Parameter(
            displayName="HUC8",
            name="huc8",
            datatype="GPString",
            parameterType="Required",
            direction="Input"
        ); p_huc8.filter.type = "ValueList"
        p_gdb   = arcpy.Parameter(
            displayName="Project Geodatabase",
            name="project_gdb",
            datatype="Workspace",           # ← was DEWorkspace, Workspace is the safest for GDBs in a .pyt.
            parameterType="Required",
            direction="Input"
        )
        p_basename  = arcpy.Parameter(
            displayName="Basename (from previous steps)",
            name="basename",
            datatype="GPString",
            parameterType="Derived",
            direction="Input"
        )
        p_outfds= arcpy.Parameter(  
            displayName="Output Feature Dataset",
            name="output_feature_dataset",
            datatype="DEFeatureDataset",
            parameterType="Derived",
            direction="Input"
        )
        p_500   = arcpy.Parameter(
            displayName="Flood_500yr",
            name="flood_500yr",
            datatype="DEFeatureClass",
            parameterType="Derived",
            direction="Input"
        )
        p_100   = arcpy.Parameter(
            displayName="Flood_100yr",
            name="flood_100yr",
            datatype="DEFeatureClass",
            parameterType="Derived",
            direction="Input"
        )
        p_10    = arcpy.Parameter(
            displayName="Flood_10yr",
            name="flood_10yr",
            datatype="DEFeatureClass",
            parameterType="Derived",
            direction="Input"
        )

        # Optional envs
        p_parallel = arcpy.Parameter(
            displayName="Parallel Processing Factor",
            name="parallel_factor",
            datatype="GPString",
            parameterType="Optional",
            direction="Input",
            category="Environment Settings"
        ); p_parallel.value = "100%"
        p_overwrite= arcpy.Parameter(
            displayName="Overwrite Outputs",
            name="overwrite_outputs",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input",
            category="Environment Settings"
        ); p_overwrite.value = True
        p_use_pairwise = arcpy.Parameter(
            displayName="Use pairwise tools where possible.\n(Faster option but may fail on complex geometries.)",
            name="use_pairwise",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input",
            category="Debug"
        ); p_use_pairwise.value = False
        p_force_buf0 = arcpy.Parameter(
            displayName="Force 0-buffer pre-clean before erase/merge.\n(Slower option but topologically safer.)",
            name="force_buffer0",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input",
            category="Debug"
        ); p_force_buf0.value = False
        p_keep_inter = arcpy.Parameter(
            displayName="Keep intermediate datasets",
            name="keep_intermediates",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input",
            category="Debug"
        ); p_keep_inter.value = False
        p_validate = arcpy.Parameter(
            displayName="Validate topology after build",
            name="validate_topology",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input"
        ); p_validate.value = False  # default OFF
        
        # Derived
        p_finalfds = arcpy.Parameter(
            displayName="Final Feature Dataset",
            name="final_feature_dataset",
            datatype="DEFeatureDataset",
            parameterType="Derived",
            direction="Output"
        )
        p_topo     = arcpy.Parameter(
            displayName="Final Topology",
            name="final_topology",
            datatype="DETopology",
            parameterType="Derived",
            direction="Output"
        )
        p_merged   = arcpy.Parameter(
            displayName="Existing_Future_<Basename>_500yr_100yr_10yr",
            name="ef_merged",
            datatype="DEFeatureClass",
            parameterType="Derived",
            direction="Output"
        )

        return [p_lut, p_watershed, p_huc8, p_gdb, p_basename, p_outfds, p_500, p_100, p_10,
                p_parallel, p_overwrite,
                p_use_pairwise, p_force_buf0, p_keep_inter, p_validate, 
                p_finalfds, p_topo, p_merged]

    def updateParameters(self, params):
        (p_lut, p_watershed, p_huc8, p_gdb, p_basename, p_outfds, p_500, p_100, p_10) = params[0:9]

        # Cascade updates based on lookup table
        if p_lut.altered and p_lut.value:
            _updateWatershedAndHUC8Filters(p_lut, p_watershed, p_huc8)

        # Update path values based on parameter values
        if p_watershed.value and p_huc8.value:
            _updateDerivedPaths(p_gdb, p_watershed, p_huc8, p_basename, p_outfds, p_500, p_100, p_10)
        return

    def execute(self, params, messages):
        (p_lut, p_watershed, p_huc8, p_gdb, p_basename, p_outfds, p_500, p_100, p_10,
         p_parallel, p_overwrite, 
         p_validate, p_use_pairwise, p_force_buf0, p_keep_inter,
         p_finalfds, p_topo, p_merged) = params
        

        apply_envs(p_parallel.valueAsText or "100%",
                   p_overwrite.value,
                   )
        # Basename: prefer explicit; else derive from WS/HUC8
        basename = p_basename.valueAsText or (compute_basename(p_watershed.valueAsText, p_huc8.valueAsText) if (p_watershed.value and p_huc8.value) else None)
        if not basename:
            raise RuntimeError("Basename is required (either supply it or choose Watershed + HUC8).")
        p_basename.value = basename

        gdb = p_gdb.valueAsText
        out_fds = os.path.join(gdb, f"Output_{basename}")
        p_500 = os.path.join(out_fds, f"Existing_Future_{basename}_500yr")
        p_100 = os.path.join(out_fds, f"Existing_Future_{basename}_100yr")
        p_10  = os.path.join(out_fds, f"Existing_Future_{basename}_10yr")

        r = combine_floods(
            gdb=p_gdb.valueAsText,
            basename=p_basename.valueAsText,
            out_fds=p_outfds.valueAsText,
            fc_500=p_500.valueAsText,
            fc_100=p_100.valueAsText,
            fc_10=p_10.valueAsText,
            validate_topology=bool(p_validate.value),
            use_pairwise=bool(p_use_pairwise.value),
            force_buffer0=bool(p_force_buf0.value),
            keep_intermediates=bool(p_keep_inter.value)
            )

        p_finalfds.value = r["final_fds"]
        p_topo.value     = r["topology"]
        p_merged.value   = r["ef_merged"]

# ---------------------------
# Tool 4: RunFloodQuiltPipeline (optimized)
# ---------------------------
class RunFloodQuiltPipeline(object):
    """
    One-click pipeline that executes:
        1) SplitFloodFreqs  -> Split by 10yr/100yr/500yr and converts to single-part
        2) MergeFloodFreqs  -> Dissolved output by 10yr/100yr/500yr
        3) CombineFloodFreqs-> Final feature dataset with merged existing and future layer and associated topology

    Improvements:
      - Centralized env application (parallel, overwrite, workspace/scratch/OCS)
      - Validates inputs and fields early; clearer messages and timings
      - Minimizes disk I/O (uses 'in_memory' for some intermediates)
      - Reuses Merge tool helper methods (grouping, SR pick, FD creation)
      - Optional skip flags for iterative debugging/reruns
    """
    def __init__(self):
        self.label = "Run Flood Quilt Pipeline"
        self.description = "End-to-end pipeline with integrated steps, robust validation, and performance-lean intermediates."
        self.canRunInBackground = False

    def getParameterInfo(self):
        # Selector
        p_lut = arcpy.Parameter(
            displayName="Lookup Table (Watershed,HUC8)",
            name="lookup_table",
            datatype="Table View",
            parameterType="Required",
            direction="Input"
        )
        ''' # Future enhancement: allow user to specify watershed and HUC fields if needed
        p_watershed_field = arcpy.Parameter( 
            displayName="Watershed Field",
            name="ws_field",
            datatype="Field",
            parameterType="Required",
            direction="Input"
        ); p_watershed_field.parameterDependencies = [p_lut.name]
        
        p_huc_field = arcpy.Parameter(
            displayName="HUC Field",
            name="huc_field",
            datatype="Field",
            parameterType="Required",
            direction="Input"
        ); p_huc_field.parameterDependencies = [p_lut.name]
        '''
        p_watershed  = arcpy.Parameter(
            displayName="Watershed",
            name="watershed",
            datatype="GPString",
            parameterType="Required",
            direction="Input"
        ); p_watershed.filter.type="ValueList"
        p_huc8 = arcpy.Parameter(
            displayName="HUC8",
            name="huc8",
            datatype="GPString",
            parameterType="Required",
            direction="Input"
        ); p_huc8.filter.type="ValueList"

        # Common
        p_gdb = arcpy.Parameter(
            displayName="Project Geodatabase",
            name="project_gdb",
            datatype="Workspace",
            parameterType="Required",
            direction="Input"
        )
        p_sr  = arcpy.Parameter(
            displayName="Target Spatial Reference",
            name="target_sr",
            datatype="GPSpatialReference", # consider GPCoordinateSystem
            parameterType="Required",
            direction="Input"
        )

        # Split inputs
        p_in_exist  = arcpy.Parameter(
            displayName="Existing Flood Feature Class",
            name="input_feature_existing",
            datatype="DEFeatureClass",
            parameterType="Required",
            direction="Input"
        )
        p_in_future = arcpy.Parameter(
            displayName="Future Flood Feature Class",
            name="input_feature_future",
            datatype="DEFeatureClass",
            parameterType="Required",
            direction="Input"
        )
        p_splitfield= arcpy.Parameter(
            displayName="Split Field (frequency)",
            name="split_field",
            datatype="GPString",
            parameterType="Optional",
            direction="Input"
        ); p_splitfield.value="FLOOD_FREQ"

        # Merge inputs
        ''' # not used, paths derived internally '''

        # Envs
        p_parallel = arcpy.Parameter(
            displayName="Parallel Processing Factor",
            name="parallel_factor",
            datatype="GPString",
            parameterType="Optional",
            direction="Input",
            category="Environment Settings"
        ); p_parallel.value = "100%"
        p_overwrite= arcpy.Parameter(
            displayName="Overwrite Outputs",
            name="overwrite_outputs",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input",
            category="Environment Settings"
        ); p_overwrite.value = True
        p_use_pairwise = arcpy.Parameter(
            displayName="Use pairwise tools where possible.\n(Faster option but may fail on complex geometries.)",
            name="use_pairwise",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input",
            category="Debug"
        ); p_use_pairwise.value = False
        p_force_buf0 = arcpy.Parameter(
            displayName="Force 0-buffer pre-clean before erase/merge.\n(Slower option but topologically safer.)",
            name="force_buffer0",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input",
            category="Debug"
        ); p_force_buf0.value = False
        p_keep_inter = arcpy.Parameter(
            displayName="Keep intermediate datasets",
            name="keep_intermediates",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input",
            category="Debug"
        ); p_keep_inter.value = False

        # Derived finals
        p_basename  = arcpy.Parameter(
            displayName="Basename (derived)",
            name="basename",
            datatype="GPString",
            parameterType="Derived",
            direction="Output"
        )
        p_finalfds  = arcpy.Parameter(
            displayName="Final Feature Dataset",
            name="final_feature_dataset",
            datatype="DEFeatureDataset",
            parameterType="Derived",
            direction="Output"
        )
        p_topology  = arcpy.Parameter(
            displayName="Final Topology",
            name="final_topology",
            datatype="DETopology",
            parameterType="Derived",
            direction="Output"
        )
        p_efmerged  = arcpy.Parameter(
            displayName="Existing_Future_<Basename>_500yr_100yr_10yr",
            name="ef_merged",
            datatype="DEFeatureClass",
            parameterType="Derived",
            direction="Output"
        )
        p_validate_topo = arcpy.Parameter(
            displayName="Validate topology after build",
            name="validate_topology",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input"
        ); p_validate_topo.value = False  # default OFF
        
        # Optional skip flags (feature enhancement)
        p_skip_split  = arcpy.Parameter(
            displayName="Skip Split step",
            name="skip_split",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input",
            category="Debug Options"
        ); p_skip_split.value = False
        p_skip_merge  = arcpy.Parameter(
            displayName="Skip Merge step",
            name="skip_merge",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input",
            category="Debug Options"
        ); p_skip_merge.value = False
        p_skip_combine= arcpy.Parameter(
            displayName="Skip Combine step",
            name="skip_combine",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input",
            category="Debug Options"
        ); p_skip_combine.value = False

        return [p_lut, p_watershed, p_huc8, p_gdb, p_sr,
                p_in_exist, p_in_future, p_splitfield,
                p_parallel, p_overwrite, 
                p_use_pairwise, p_force_buf0, p_keep_inter, 
                p_basename, p_finalfds, p_topology, p_efmerged,
                p_validate_topo,
                p_skip_split, p_skip_merge, p_skip_combine]

    def updateParameters(self, params):
        (p_lut, p_watershed, p_huc8, p_gdb) = params[0:3]
        p_basename = params[13] # derived

        # Cascade updates based on lookup table
        if p_lut.altered and p_lut.value:
            self._updateWatershedAndHUC8Filters(p_lut, p_watershed, p_huc8)
    
    ''' # Future enhancement: allow user to specify watershed and HUC fields if needed'''

    # ---- Internal helpers for this module's execution ----
    def _validate_field(self, fc, field_name):
        field_names = [f.name for f in arcpy.ListFields(fc)]
        if field_name not in field_names:
            raise RuntimeError(f"Required field '{field_name}' not found in {fc}.")

    def _ensure_fd(self, gdb, name, sr):
        out = os.path.join(gdb, name)
        if not arcpy.Exists(out):
            arcpy.management.CreateFeatureDataset(gdb, name, sr)
        return out

    def _repair(self, fc):
        return arcpy.management.RepairGeometry(fc, validation_method="OGC").getOutput(0)

    def execute(self, params, messages):
        (p_lut, p_ws, p_huc8, p_gdb, p_sr,
        p_in_exist, p_in_future, p_splitfield,
        p_table, p_field,
        p_parallel, p_overwrite, 
        p_use_pairwise, p_force_buf0, p_keep_inter,
        p_basename, p_finalfds, p_topology, p_efmerged,
        p_validate_topo,
        p_skip_split, p_skip_merge, p_skip_combine) = params

        apply_envs(p_parallel.valueAsText or "100%",
                p_overwrite.value,
                )

        log = RunLog()
        gdb     = p_gdb.valueAsText
        sr_obj  = p_sr.value
        basename = compute_basename(p_ws.valueAsText, p_huc8.valueAsText) if (p_ws.value and p_huc8.value) else p_basename.valueAsText
        if not basename: raise RuntimeError("Basename is required (choose Watershed+HUC8 or provide Basename).")
        p_basename.value = basename

        # -------- 1) SPLIT (shared helper) --------
        if not p_skip_split.value:
            _msg("1) SPLIT SUBMODULE: Creating existing and future datasets and splitting by frequency…")
            out = split_floods(gdb, basename, sr_obj, p_in_exist.valueAsText, p_in_future.valueAsText,
                            split_field=(p_splitfield.valueAsText or "FLOOD_FREQ"))
            log.add("Split Existing FDS", out["fds_exist"])
            log.add("Split Future FDS",   out["fds_future"])
        else:
            arcpy.AddMessage("SPLIT SUBMODULE SKIPPED BY USER.", warn=True)
            out = {}

        # -------- 2) MERGE (shared helper) --------
        if not p_skip_merge.value:
            _msg("2) MERGE SUBMODULE: Merging and dissolving by FLOOD_FREQ into output dataset.")
            paths = read_fc_paths_from_table(p_table.valueAsText, p_field.valueAsText or "FCPath")
            if len(paths) < 3: raise RuntimeError("Merge step: expected at least 3 rows with FC paths.")
            buckets = group_by_suffix(paths)
            sr = choose_spatial_ref_from_inputs(paths, sr_obj)
            rmerge = merge_floods(project_gdb=gdb,
                                basename=basename,
                                buckets=buckets,
                                target_sr=sr,
                                use_pairwise=bool(p_use_pairwise.value),
                                force_buffer0=bool(p_force_buf0.value),
                                keep_intermediates=bool(p_keep_inter.value))
            out_fds = rmerge["out_fds"]
            out_500 = rmerge["out_500"]; out_100 = rmerge["out_100"]; out_010 = rmerge["out_010"]
        else:
            arcpy.AddMessage("MERGE SUBMODULE SKIPPED BY USER.", warn=True)
            out_fds = os.path.join(gdb, f"Output_{basename}")
            out_500 = os.path.join(out_fds, f"Existing_Future_{basename}_500yr")
            out_100 = os.path.join(out_fds, f"Existing_Future_{basename}_100yr")
            out_010 = os.path.join(out_fds, f"Existing_Future_{basename}_10yr")

        
            if not any(arcpy.Exists(fc) for fc in (out_500, out_100, out_010)):
                raise RuntimeError("Merge skipped but no existing Flood_* layers found in Output_<Basename>.")

        # -------- 3) COMBINE (shared helper) --------
        if not p_skip_combine.value:
            _msg("3) COMBINE SUBMODULE: Repairing, erasing overlaps, merging, and building topology for final features.")
            rcomb = combine_floods(
                gdb=gdb, basename=basename, out_fds=out_fds,
                fc_500=out_500, fc_100=out_100, fc_10=out_010,
                validate_topology=bool(p_validate_topo.value),
                use_pairwise=bool(p_use_pairwise.value),
                force_buffer0=bool(p_force_buf0.value),
                keep_intermediates=bool(p_keep_inter.value)
            )
            p_finalfds.value = rcomb["final_fds"]
            p_topology.value = rcomb["topology"]
            p_efmerged.value = rcomb["ef_merged"]
            log.add("Final FDS", rcomb["final_fds"])
            log.add("Final Topology", rcomb["topology"])
            log.add("EF merged", rcomb["ef_merged"])
        else:
            arcpy.AddMessage("COMBINE SUBMODULE SKIPPED BY USER.", warn=True)

        # Summary
        log.dump("Pipeline summary")
