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
- 2025-09-26: Initial version, based on FloodQuiltSplitModule.py, FloodQuiltMergeModule.py, and FloodQuiltCombineModule.py
- 2025-09-29: Added Pipeline tool to run all steps in sequence. Enabled topology validation as an option.
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
from collections import defaultdict

# ---------------------------
# Global/Shared helpers
# ---------------------------

def _msg(msg, warn=False): 
    """Helper to print a message or warning."""
    (arcpy.AddWarning if warn else arcpy.AddMessage)(msg)

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

def split_floods(gdb: str, basename: str, spatial_ref, existing_fc: str, future_fc: str, split_field: str = "FLOOD_FREQ"):
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
    _msg(f"Splitting existing feature class by {split_field}…")
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
    _msg("Converting multipart to singlepart for existing outputs…")
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
    _msg(f"Splitting future feature class by {split_field}…")
    arcpy.analysis.SplitByAttributes(future_fc, fds_future, [split_field])

    F10 = _pick(fds_future, ["T10", "10yr", "10", "YR10", "yr10"])
    F1 = _pick(fds_future, ["T1", "100yr", "100", "YR100", "yr100"])
    F5 = _pick(fds_future, ["T0_2", "500yr", "500", "YR500", "yr500"])

    # Canonical single-part outputs for future
    future_10 = os.path.join(fds_future, f"Future_{basename}_10yr") if F10 else ""
    future_100 = os.path.join(fds_future, f"Future_{basename}_100yr") if F1 else ""
    future_500 = os.path.join(fds_future, f"Future_{basename}_500yr") if F5 else ""

    # Multipart → Singlepart for future
    _msg("Converting multipart to singlepart for future outputs…")
    for src, dst in ((F10, future_10), (F1, future_100), (F5, future_500)):
        if src and dst:
            arcpy.management.MultipartToSinglepart(src, dst)

    # Delete raw split outputs if converted for future
    to_delete_future = [p for p in (F10, F1, F5) if p and arcpy.Exists(p)]
    if to_delete_future:
        arcpy.management.Delete(to_delete_future)
    _msg("Future feature class processing complete.")

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

def combine_floods( 
    gdb: str,
    basename: str,
    out_feature_dataset: str,      # usually Output_<Basename>
    fc_500: str,
    fc_100: str,
    fc_10: str,
    out_cs=None,
    validate_topology: bool = False   # toggle validation
):
    """
    Core Combine step used by Tool 3 (CombineFloodFreqs) and Tool 4 (Pipeline).

    Steps:
      - Create Final_<Basename> feature dataset (SR from 10yr/100yr/500yr or env/output CS)
      - Repair geometry (OGC)
      - Erase: 500yr-100yr  → <basename>_500yr_Erased (in Output_<Basename>)
               100yr-10yr   → <basename>_100yr_Erased (in Output_<Basename>)
      - Merge: [10yr_repaired, 500yr_erased, 100yr_erased]
               → Final_<Basename>/Existing_Future_<Basename>_500yr_100yr_10yr
      - Create topology + optional ValidateTopology
    """
    import arcpy, os

    final_fds_name = f"Final_{basename}"
    final_fds = os.path.join(gdb, final_fds_name)
    _msg(f"Combining into feature dataset: {final_fds}")

    # Choose SR for Final FDS
    sr_final = None
    for fc in (fc_10, fc_100, fc_500):
        if fc and arcpy.Exists(fc):
            dsr = arcpy.Describe(fc).spatialReference
            if dsr and dsr.factoryCode:
                sr_final = dsr
                break
    if sr_final is None:
        sr_final = out_cs or arcpy.env.outputCoordinateSystem
    _msg(f"Using spatial reference: {sr_final.name}")

    if not arcpy.Exists(final_fds):
        arcpy.management.CreateFeatureDataset(gdb, final_fds_name, sr_final)
    _msg("Final feature dataset created")

    # Repair
    _msg("Repairing geometries…")
    rep_500 = arcpy.management.RepairGeometry(fc_500, validation_method="OGC").getOutput(0) if (fc_500 and arcpy.Exists(fc_500)) else None
    rep_100 = arcpy.management.RepairGeometry(fc_100, validation_method="OGC").getOutput(0) if (fc_100 and arcpy.Exists(fc_100)) else None
    rep_010 = arcpy.management.RepairGeometry(fc_10,  validation_method="OGC").getOutput(0) if (fc_10  and arcpy.Exists(fc_10))  else None

    # Erase overlaps, stored in Output_<Basename>
    _msg("Erasing overlaps…")
    erased_500 = os.path.join(out_feature_dataset, f"{basename}_500yr_Erased")
    erased_100 = os.path.join(out_feature_dataset, f"{basename}_100yr_Erased")
    erased_10 = os.path.join(out_feature_dataset, f"{basename}_10yr_Erased")

    _msg("Erasing 500yr-100yr…")
    if rep_500 and rep_100:
        arcpy.analysis.Erase(rep_500, rep_100, erased_500)
    else:
        erased_500 = ""
    _msg("Erasing 100yr-10yr…")
    if rep_100 and rep_010:
        arcpy.analysis.Erase(rep_100, rep_010, erased_100)
    else:
        erased_100 = ""

    sp_500 = os.path.join(out_feature_dataset, f"{basename}_500yr_Erased_Singlepart")
    sp_100 = os.path.join(out_feature_dataset, f"{basename}_100yr_Erased_Singlepart")
    sp_010 = os.path.join(out_feature_dataset, f"{basename}_10yr_Erased_Singlepart")

    # Multipart → Singlepart for erased outputs
    for src, dst in ((erased_500, sp_500), (erased_100, sp_100), (rep_010, sp_010)):
        if src and dst and arcpy.Exists(src):
            arcpy.management.MultipartToSinglepart(src, dst)

    # Use singlepart outputs if they exist
    inputs = [fc for fc in (sp_010 if arcpy.Exists(sp_010) else rep_010,
                            sp_100 if arcpy.Exists(sp_100) else erased_100,
                            sp_500 if arcpy.Exists(sp_500) else erased_500)
              if fc and arcpy.Exists(fc)]
    if not inputs:
        raise RuntimeError("Combine step: nothing to merge. Ensure Merge outputs exist.")
    
    # Merge thematic stack
    _msg("Merging layers into final…")
    ef_merged = os.path.join(final_fds, f"Existing_Future_{basename}_500yr_100yr_10yr")
    arcpy.management.Merge(inputs, ef_merged)
    _msg(f"Merge sucessful: {ef_merged}")

    
    '''# Merge thematic stack
    _msg("Merging layers into final…")
    ef_merged = os.path.join(final_fds, f"Existing_Future_{basename}_500yr_100yr_10yr")
    inputs = [fc for fc in (rep_010, erased_100 if erased_100 else None, erased_500 if erased_500 else None,) if fc and arcpy.Exists(fc)]
    if not inputs:
        raise RuntimeError("Combine step: nothing to merge. Ensure Merge outputs exist.")
    arcpy.management.Merge(inputs, ef_merged)
    _msg(f"Merge sucessful: {ef_merged}")'''

    arcpy.management.RepairGeometry(ef_merged, validation_method="OGC")
    
    # Topology
    _msg("Creating topology…")
    topo = arcpy.management.CreateTopology(final_fds, f"Final_{basename}_Topology").getOutput(0)
    arcpy.management.AddFeatureClassToTopology(topo, ef_merged, 1, 1)
    arcpy.management.AddRuleToTopology(topo, "Must Not Overlap (Area)", ef_merged)
    _msg("Topology creation sucessful.")

    # Optional validation (can be slow)
    if validate_topology:
        arcpy.AddMessage("Validating topology (Full Extent)…")
        arcpy.management.ValidateTopology(topo, "Full_Extent")

    return {
        "final_fds": final_fds,
        "topology": topo,
        "ef_merged": ef_merged,
        "erased_500": erased_500,
        "erased_100": erased_100
    }

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
        p_ws  = arcpy.Parameter(
            displayName="Watershed",
            name="watershed",
            datatype="GPString",
            parameterType="Required",
            direction="Input"
        ); p_ws.filter.type = "ValueList"
        p_h8  = arcpy.Parameter(
            displayName="HUC8",
            name="huc8",
            datatype="GPString",
            parameterType="Required",
            direction="Input"
        ); p_h8.filter.type = "ValueList"

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
        '''
        p_wsenv    = arcpy.Parameter(       # not used, envs applied globally
            displayName="Workspace",
            name="env_workspace",
            datatype="Workspace",           # ← was DEWorkspace, Workspace is the safest for GDBs in a .pyt.
            parameterType="Optional",
            direction="Input",
            category="Environment Settings"
        )
        p_scratch  = arcpy.Parameter(       # not used, envs applied globally
            displayName="Scratch Workspace",
            name="env_scratch",
            datatype="Workspace",           # ← was DEWorkspace, Workspace is the safest for GDBs in a .pyt.
            parameterType="Optional",
            direction="Input",
            category="Environment Settings"
        )
        p_ocs      = arcpy.Parameter(       # not used, envs applied globally
            displayName="Output Coordinate System",
            name="env_ocs",
            datatype="GPSpatialReference",
            parameterType="Optional",
            direction="Input",
            category="Environment Settings"
        )
        '''
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

        return [p_lut, p_ws, p_h8, p_gdb, p_sr, p_in_exist, p_in_future, p_splitfield,
            p_parallel, p_overwrite, #p_wsenv, p_scratch, p_ocs, # not used, envs applied globally
            p_basename, p_fds_exist, p_exist_10, p_exist_100, p_exist_500,
            p_fds_future, p_future_10, p_future_100, p_future_500]

    def updateParameters(self, params):
        (p_lut, p_ws, p_h8, *_rest) = params[:3+7]  # just to keep readable
        p_basename = params[10] # derived

        # Cascade
        if p_lut.altered and p_lut.value:
            lookup = build_ws_huc_lookup(p_lut.valueAsText)
            p_ws.filter.list = sorted(lookup.keys())
            if p_ws.value and p_ws.value in lookup:
                p_h8.filter.list = lookup[p_ws.value]
            else:
                p_h8.filter.list = []

        if p_ws.value and p_h8.value:
            p_basename.value = compute_basename(p_ws.value, p_h8.value)
        return

    def execute(self, params, messages):
        (p_lut, p_ws, p_h8, p_gdb, p_sr, p_in_exist, p_in_future, p_splitfield,
         p_parallel, p_overwrite, #p_wsenv, p_scratch, p_ocs, # not used, envs applied globally
         p_basename, p_fds_exist, p_exist_10, p_exist_100, p_exist_500,
         p_fds_future, p_future_10, p_future_100, p_future_500) = params

        apply_envs(p_parallel.valueAsText or "100%",
                   p_overwrite.value,  # currently not overwriting instead appending 0 and end of name
                   #p_wsenv.valueAsText,
                   #p_scratch.valueAsText,
                   #p_ocs.value
                   )

        gdb         = p_gdb.valueAsText
        spatial_ref = p_sr.value
        in_exist    = p_in_exist.valueAsText
        in_future   = p_in_future.valueAsText
        split_field = (p_splitfield.valueAsText or "FLOOD_FREQ")

        basename = compute_basename(p_ws.valueAsText, p_h8.valueAsText)
        arcpy.SetParameterAsText(10, basename)  # derived
        _msg(f"Basename: {basename}")
        _msg(f"Splitting by field: {split_field}")

        # Shared core
        out = split_floods(gdb, 
                           basename, 
                           spatial_ref, 
                           in_exist, 
                           in_future, 
                           split_field
                           )
        _msg(f"Split complete. Outputs:\n{out}")   

        # Push derived outputs
        p_fds_exist.value  = out["fds_exist"]
        p_fds_future.value = out["fds_future"]
        p_exist_10.value   = out["exist_10"]
        p_exist_100.value  = out["exist_100"]
        p_exist_500.value  = out["exist_500"]
        p_future_10.value  = out["future_10"]
        p_future_100.value = out["future_100"]
        p_future_500.value = out["future_500"]

        # Paths to feature classes that were modified
        fc_paths = {
            "Existing_10yr": out["exist_10"],
            "Existing_100yr": out["exist_100"],
            "Existing_500yr": out["exist_500"],
            "Future_10yr": out["future_10"],
            "Future_100yr": out["future_100"],
            "Future_500yr": out["future_500"],
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

        _msg("Split Flood Frequencies completed.")

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
        p_ws  = arcpy.Parameter(
            displayName="Watershed",
            name="watershed",
            datatype="GPString",
            parameterType="Required",
            direction="Input"
        ); p_ws.filter.type = "ValueList"
        p_h8  = arcpy.Parameter(
            displayName="HUC8",
            name="huc8",
            datatype="GPString",
            parameterType="Required",
            direction="Input"
        ); p_h8.filter.type = "ValueList"
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
        p_base  = arcpy.Parameter(         # Check if being passed correctly
            displayName="Basename (from Split tool)",
            name="basename",
            datatype="GPString",
            parameterType="Required",
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
        '''
        p_wsenv    = arcpy.Parameter(       # not used, envs applied globally
            displayName="Workspace",
            name="env_workspace",
            datatype="Workspace",           # ← was DEWorkspace, Workspace is the safest for GDBs in a .pyt.
            parameterType="Optional",
            direction="Input",
            category="Environment Settings"
        )
        p_scratch  = arcpy.Parameter(       # not used, envs applied globally
            displayName="Scratch Workspace",
            name="env_scratch",
            datatype="Workspace",           # ← was DEWorkspace, Workspace is the safest for GDBs in a .pyt.
            parameterType="Optional",
            direction="Input",
            category="Environment Settings"
        )
        p_ocs      = arcpy.Parameter(       # not used, envs applied globally
            displayName="Output Coordinate System",
            name="env_ocs",
            datatype="GPSpatialReference",
            parameterType="Optional",
            direction="Input",
            category="Environment Settings"
        )
        '''

        # Derived
        p_outfds = arcpy.Parameter(
            displayName="Output Feature Dataset",
            name="out_feature_dataset",
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

        return [p_lut, p_ws, p_h8, p_table, p_field, p_gdb, p_base, p_sr,
                p_parallel, p_overwrite, #p_wsenv, p_scratch, p_ocs,    # not used, envs applied globally
                p_outfds, p_out500, p_out100, p_out10]

    def updateParameters(self, params):
        p_lut, p_ws, p_h8 = params[0], params[1], params[2]  # just to keep readable
        p_basename = params[6] # derived

        # Cascade
        if p_lut.altered and p_lut.value:
            lookup = build_ws_huc_lookup(p_lut.valueAsText)
            p_ws.filter.list = sorted(lookup.keys())
            if p_ws.value and p_ws.value in lookup:
                p_h8.filter.list = lookup[p_ws.value]
            else:
                p_h8.filter.list = []

        if p_ws.value and p_h8.value:
            p_basename.value = compute_basename(p_ws.value, p_h8.value)
        return

    def _scratch_gdb(self):
        sgdb = arcpy.env.scratchGDB
        if sgdb and arcpy.Exists(sgdb):
            return sgdb
        # Fallback: create a temp file GDB in the scratchFolder (or Workspace)
        base_folder = (arcpy.env.scratchFolder or arcpy.env.workspace or arcpy.env.scratchWorkspace or os.getcwd())
        name = arcpy.CreateUniqueName("fgdb_", base_folder)  # generates unique folder-like name
        if not name.lower().endswith(".gdb"):
            name = f"{name}.gdb"
        arcpy.management.CreateFileGDB(os.path.dirname(name), os.path.basename(name))
        return name

    def _same_sr(self, fc, sr):
        try:
            dsr = arcpy.Describe(fc).spatialReference
            return (dsr and sr and dsr.factoryCode == sr.factoryCode and dsr.name == sr.name)
        except Exception:
            return False

    def _normalize_inputs(self, fc_list, target_sr):
        """Project + repair each input into scratch GDB to enforce common SR/XY properties."""
        out = []
        scratch = self._scratch_gdb()
        for i, fc in enumerate(fc_list, 1):
            if not arcpy.Exists(fc):
                continue
            tmp = arcpy.CreateUniqueName(f"norm_{i}", scratch)
            if self._same_sr(fc, target_sr):
                # Copy as-is to ensure consistent storage + tolerance
                arcpy.management.CopyFeatures(fc, tmp)
            else:
                arcpy.management.Project(fc, tmp, target_sr)
            # Repair (delete null shapes if present)
            try:
                arcpy.management.RepairGeometry(tmp, validation_method="OGC", delete_null="DELETE_NULL")
            except TypeError:
                # older Pro builds don’t have delete_null; fall back
                arcpy.management.RepairGeometry(tmp, validation_method="OGC")
            out.append(tmp)
        return out

    def _preclean_fc(self, in_fc):
        """Make a cleaned copy in scratch GDB to stabilize topology before dissolve."""
        scratch = self._scratch_gdb()
        base = os.path.basename(in_fc)
        clean = arcpy.CreateUniqueName(f"clean_{base}", scratch)
        arcpy.management.CopyFeatures(in_fc, clean)
        try:
            arcpy.management.RepairGeometry(clean, validation_method="OGC", delete_null="DELETE_NULL")
        except TypeError:
            arcpy.management.RepairGeometry(clean, validation_method="OGC")
        # Optional polygon 0-buffer to fix self-intersections; no-op for lines
        gtype = arcpy.Describe(clean).shapeType
        if gtype.lower() == "polygon":
            buf = arcpy.CreateUniqueName(f"buf0_{base}", scratch)
            arcpy.analysis.Buffer(clean, buf, "0 Feet", dissolve_option="NONE")
            clean = buf
        return clean

    def _safe_dissolve(self, in_fc, out_fc, freq_field):
        """Dissolve with a fallback path if 160196 occurs."""
        # Choose line/area option properly
        gtype = arcpy.Describe(in_fc).shapeType
        kwargs = dict(in_features=in_fc, out_feature_class=out_fc, dissolve_field=freq_field,
                    statistics_fields=[], multi_part="MULTI_PART")
        if gtype.lower() == "polyline":
            kwargs["unsplit_lines"] = "DISSOLVE_LINES"
        try:
            arcpy.management.Dissolve(**kwargs)
            return out_fc
        except arcpy.ExecuteError as e:
            msg = arcpy.GetMessages()
            if "160196" in msg or "Invalid Topology" in msg:
                arcpy.AddWarning("Dissolve failed with ERROR 160196; attempting pre-clean + retry…")
                clean = self._preclean_fc(in_fc)
                kwargs["in_features"] = clean
                arcpy.management.Dissolve(**kwargs)
                return out_fc
            raise
    
    def read_fc_paths_from_table(self, table, path_field): # returns list of strings
        paths = []
        with arcpy.da.SearchCursor(table, [path_field]) as cur:
            for (p,) in cur:
                if p:
                    paths.append(str(p))
        return paths

    def group_by_suffix(self, paths): # returns dict: suffix -> list of paths
        buckets = defaultdict(list)
        for p in paths:
            base = os.path.splitext(os.path.basename(p))[0]
            m = self.SUFFIX_RE.search(base)
            if m:
                key = m.group(1).lower()
                buckets[key].append(p)
            else:
                _msg(f"Skipped (no suffix match): {p}", warn=True)
        return buckets

    def pick_spatial_ref(self, paths, fallback_sr=None): # returns arcpy.SpatialReference
        return choose_spatial_ref_from_inputs(paths, fallback_sr)

    def create_feature_dataset(self, gdb, basename, sr): # returns path to FDS
        fds_name = f"Output_{basename}"
        out_fds = os.path.join(gdb, fds_name)
        if not arcpy.Exists(out_fds):
            _msg(f"Creating feature dataset: {out_fds}")
            arcpy.management.CreateFeatureDataset(gdb, fds_name, sr)
        else:
            _msg(f"Feature Dataset exists, using: {out_fds}")
        return out_fds

    def ensure_field_exists(self, in_fc, field_name): # raises RuntimeError if not found
        if field_name not in [f.name for f in arcpy.ListFields(in_fc)]:
            raise RuntimeError(f"Field '{field_name}' not found in {in_fc}.")

    def merge_group(self, fc_list, out_fc, target_sr=None):
        if len(fc_list) == 0:
            raise RuntimeError(f"No inputs for {out_fc}")
        # Normalize → Merge
        norm = self._normalize_inputs(fc_list, target_sr)
        arcpy.AddMessage(f"Merging {len(norm)} normalized inputs -> {out_fc}")
        arcpy.management.Merge(norm, out_fc)
        # Pre-clean merged temp to stabilize topology for dissolve later
        return out_fc
    
    def dissolve_on_freq(self, in_fc, out_fc, freq_field=DISSOLVE_FIELD): # returns path to dissolved FC
        self.ensure_field_exists(in_fc, freq_field)
        _msg(f"Dissolving {in_fc} by {freq_field} -> {out_fc}")
        self._safe_dissolve(in_fc, out_fc, freq_field)
        return out_fc

    def execute(self, params, messages):
        (p_lut, p_ws, p_h8, p_table, p_field, p_gdb, p_base, p_sr,
        p_parallel, p_overwrite,
        p_outfds, p_out500, p_out100, p_out10) = params

        apply_envs(p_parallel.valueAsText or "100%", p_overwrite.value)

        input_table = p_table.valueAsText
        path_field  = p_field.valueAsText or "FCPath"
        project_gdb = p_gdb.valueAsText
        sr_param    = p_sr.value

        # Basename: prefer explicit input; else derive from WS/HUC8
        basename = p_base.valueAsText or (compute_basename(p_ws.valueAsText, p_h8.valueAsText) if (p_ws.value and p_h8.value) else None)
        if not basename:
            raise RuntimeError("Basename is required (either connect from Split or choose Watershed + HUC8).")
        p_base.value = basename  # keep dialog synced

        # Read & group
        paths = self.read_fc_paths_from_table(input_table, path_field)
        if len(paths) < 3:
            raise RuntimeError("Expected at least 3 rows with FC paths.")
        buckets = self.group_by_suffix(paths)
        for key in self.FREQ_KEYS:
            if not buckets.get(key):
                _msg(f"No inputs found for '{key}' group.", warn=True)

        # SR + FDS
        sr = self.pick_spatial_ref(paths, sr_param)
        out_fds = self.create_feature_dataset(project_gdb, basename, sr)

        # Merge temps
        merged_500 = os.path.join(out_fds, f"Merged_{basename}_500yr")
        merged_100 = os.path.join(out_fds, f"Merged_{basename}_100yr")
        merged_010 = os.path.join(out_fds, f"Merged_{basename}_10yr")

        # Final
        out_500 = os.path.join(out_fds, f"Existing_Future_{basename}_500yr")
        out_100 = os.path.join(out_fds, f"Existing_Future_{basename}_100yr")
        out_010 = os.path.join(out_fds, f"Existing_Future_{basename}_10yr")

        '''
        # Merge and dissolve with normalization & safe fallback
        if buckets.get("500yr"):
            self.merge_group(buckets["500yr"], merged_500, target_sr=sr)
            self.dissolve_on_freq(merged_500, out_500)
        if buckets.get("100yr"):
            self.merge_group(buckets["100yr"], merged_100, target_sr=sr)
            self.dissolve_on_freq(merged_100, out_100)
        if buckets.get("10yr"):
            self.merge_group(buckets["10yr"], merged_010, target_sr=sr)
            self.dissolve_on_freq(merged_010, out_010)
        '''

        # Merge and Dissolve (simple, without normalization or safe fallback)
        _msg("Merging and dissolving feature classes by frequency group...")
        _msg("Merging 500yr group...")
        if buckets.get("500yr"):
            arcpy.management.Merge(buckets["500yr"] + buckets.get("100yr", []) + buckets.get("10yr", []), merged_500); _msg("500yr merge complete.")
        _msg("Repairing geometry on merged 500yr...")
        if arcpy.Exists(merged_500):
            try:
                arcpy.management.RepairGeometry(merged_500, validation_method="OGC", delete_null="DELETE_NULL"); _msg("500yr repair complete.")
            except TypeError:
                arcpy.management.RepairGeometry(merged_500, validation_method="OGC"); _msg("500yr repair complete (without delete_null).")
        if arcpy.Exists(merged_500):
            try:
                # Add a new field for flood frequency if it doesn't exist
                fields = [field.name for field in arcpy.ListFields(merged_500)]
                if self.DISSOLVE_FIELD not in fields:
                    arcpy.management.AddField(merged_500, self.DISSOLVE_FIELD, "Text")
                
                # Calculate the flood frequency value as 0.2
                arcpy.management.CalculateField(merged_500, self.DISSOLVE_FIELD, "0.2", "PYTHON3")
                _msg("Flood frequency calculation complete.")
            except Exception as e:
                _msg(f"Error during flood frequency calculation: {e}")
        _msg("Dissolving 500yr group...")
        if arcpy.Exists(merged_500):
            arcpy.management.Dissolve(merged_500, out_500, dissolve_field=self.DISSOLVE_FIELD,
                                      statistics_fields=None, multi_part="MULTI_PART", unsplit_lines="DISSOLVE_LINES"); _msg("500yr dissolve complete.")
        _msg(f"500yr Feature Class: {out_500}")
        
        _msg("Merging 100yr group...")
        if buckets.get("100yr"):
            arcpy.management.Merge(buckets["100yr"] + buckets.get("10yr", []), merged_100); _msg("100yr merge complete.")
        _msg("Repairing geometry on merged 100yr...")
        if arcpy.Exists(merged_100):
            try:
                arcpy.management.RepairGeometry(merged_100, validation_method="OGC", delete_null="DELETE_NULL"); _msg("100yr repair complete.")
            except TypeError:
                arcpy.management.RepairGeometry(merged_100, validation_method="OGC"); _msg("100yr repair complete (without delete_null).")
        if arcpy.Exists(merged_100):
            try:
                # Add a new field for flood frequency if it doesn't exist
                fields = [field.name for field in arcpy.ListFields(merged_100)]
                if self.DISSOLVE_FIELD not in fields:
                    arcpy.management.AddField(merged_100, self.DISSOLVE_FIELD, "Text")
                
                # Calculate the flood frequency value as 1
                arcpy.management.CalculateField(merged_100, self.DISSOLVE_FIELD, "1", "PYTHON3")
                _msg("Flood frequency calculation complete.")
            except Exception as e:
                _msg(f"Error during flood frequency calculation: {e}")
        _msg("Dissolving 100yr group...")
        if arcpy.Exists(merged_100):
            arcpy.management.Dissolve(merged_100, out_100, dissolve_field=self.DISSOLVE_FIELD,
                                      statistics_fields=None, multi_part="MULTI_PART", unsplit_lines="DISSOLVE_LINES"); _msg("100yr dissolve complete.")
        _msg(f"100yr Feature Class: {out_100}")
        
        _msg("Merging 10yr group...")
        if buckets.get("10yr"):
            arcpy.management.Merge(buckets["10yr"], merged_010); _msg("10yr merge complete.")
        _msg("Repairing geometry on merged 10yr...")
        if arcpy.Exists(merged_010):
            try:
                arcpy.management.RepairGeometry(merged_010, validation_method="OGC", delete_null="DELETE_NULL"); _msg("10yr repair complete.")
            except TypeError:
                arcpy.management.RepairGeometry(merged_010, validation_method="OGC"); _msg("10yr repair complete (without delete_null).")
        if arcpy.Exists(merged_010):
            try:
                # Add a new field for flood frequency if it doesn't exist
                fields = [field.name for field in arcpy.ListFields(merged_010)]
                if self.DISSOLVE_FIELD not in fields:
                    arcpy.management.AddField(merged_010, self.DISSOLVE_FIELD, "Text")
                
                # Calculate the flood frequency value as 10
                arcpy.management.CalculateField(merged_010, self.DISSOLVE_FIELD, "10", "PYTHON3")
                _msg("Flood frequency calculation complete.")
            except Exception as e:
                _msg(f"Error during flood frequency calculation: {e}")
        _msg("Dissolving 10yr group...")
        if arcpy.Exists(merged_010):
            arcpy.management.Dissolve(merged_010, out_010, dissolve_field=self.DISSOLVE_FIELD,
                                      statistics_fields=None, multi_part="MULTI_PART", unsplit_lines="DISSOLVE_LINES"); _msg("10yr dissolve complete.")
        _msg(f"10yr Feature Class: {out_010}")

        '''# Cleanup temps
        for tmp in (merged_500, merged_100, merged_010): # Leaving temps for now for debugging
            if arcpy.Exists(tmp):
                arcpy.management.Delete(tmp)'''

        # Derived outputs
        p_outfds.value = out_fds
        p_out500.value = out_500 if arcpy.Exists(out_500) else ""
        p_out100.value = out_100 if arcpy.Exists(out_100) else ""
        p_out10.value  = out_010 if arcpy.Exists(out_010) else ""

        _msg("Merge & Dissolve complete.")

    """
    def execute(self, params, messages):
        (p_lut, p_ws, p_h8, p_table, p_field, p_gdb, p_base, p_sr,
         p_parallel, p_overwrite, #p_wsenv, p_scratch, p_ocs,   # not used, envs applied globally
         p_outfds, p_out500, p_out100, p_out10) = params

        apply_envs(p_parallel.valueAsText or "100%",
                   p_overwrite.value,
                   #p_wsenv.valueAsText,        # not used, envs applied globally 
                   #p_scratch.valueAsText,      # not used, envs applied globally
                   #p_ocs.value                 # not used, envs applied globally
                   )

        input_table = p_table.valueAsText
        path_field  = p_field.valueAsText or "FCPath"
        project_gdb = p_gdb.valueAsText
        basename    = p_base.valueAsText
        sr_param    = p_sr.value

        basename = p_base.valueAsText or (compute_basename(p_ws.valueAsText, p_h8.valueAsText) if (p_ws.value and p_h8.value) else None)
        if not basename:
            raise RuntimeError("Basename is required (supply it or pick Watershed + HUC8).")
        p_base.value = basename  # keep the dialog in sync

        # Read paths from table
        paths = self.read_fc_paths_from_table(input_table, path_field)
        if len(paths) < 3:
            raise RuntimeError("Expected at least 3 rows with FC paths.")

        # Group by suffix
        buckets = self.group_by_suffix(paths)
        for key in self.FREQ_KEYS:
            if not buckets.get(key):
                _msg(f"No inputs found for '{key}' group.", warn=True)

        sr = self.pick_spatial_ref(paths, sr_param)
        out_fds = self.create_feature_dataset(project_gdb, basename, sr)

        # Temporary merged paths
        merged_500 = os.path.join(out_fds, "Merged_500yr_tmp")
        merged_100 = os.path.join(out_fds, "Merged_100yr_tmp")
        merged_010 = os.path.join(out_fds, "Merged_10yr_tmp")

        # Final paths
        _msg(f"Output Feature Dataset: {out_fds}")
        out_500 = os.path.join(out_fds, f"Existing_Future_{basename}_500yr")
        out_100 = os.path.join(out_fds, f"Existing_Future_{basename}_100yr")
        out_010 = os.path.join(out_fds, f"Existing_Future_{basename}_10yr")
        
        # Merge by group
        _msg("Merging and dissolving feature classes by frequency group...")
        
        _msg("Merging 500yr group...")
        if buckets.get("500yr"):
            self.merge_group(buckets["500yr"], merged_500, target_sr=sr); _msg("500yr merge complete.")
        _msg("Dissolving 500yr group...")
        if arcpy.Exists(merged_500): 
            self.dissolve_on_freq(merged_500, out_500); _msg("500yr dissolve complete.")
            '''arcpy.management.Dissolve(merged_500, out_500, dissolve_field=self.DISSOLVE_FIELD,
                                      statistics_fields=None, multi_part="MULTI_PART", unsplit_lines="DISSOLVE_LINES"); _msg("500yr dissolve complete.")'''
        _msg(f"500yr Feature Class: {out_500}")

        _msg("Merging 100yr group...")
        if buckets.get("100yr"):
            self.merge_group(buckets["100yr"], merged_100, target_sr=sr); _msg("100yr merge complete.")
        _msg("Dissolving 100yr group...")
        if arcpy.Exists(merged_100): 
            self.dissolve_on_freq(merged_100, out_100) ;_msg("100yr dissolve complete.")
            '''arcpy.management.Dissolve(merged_100, out_100, dissolve_field=self.DISSOLVE_FIELD,
                                      statistics_fields=None, multi_part="MULTI_PART", unsplit_lines="DISSOLVE_LINES"); _msg("100yr dissolve complete.")'''
        _msg(f"100yr Feature Class: {out_100}")

        _msg("Merging 10yr group...")
        if buckets.get("10yr"):
            self.merge_group(buckets["10yr"], merged_010, target_sr=sr); _msg("10yr merge complete.")
        _msg("Dissolving 10yr group...")
        if arcpy.Exists(merged_010): 
            self.dissolve_on_freq(merged_010, out_010); _msg("10yr dissolve complete.")
            '''arcpy.management.Dissolve(merged_010, out_010, dissolve_field=self.DISSOLVE_FIELD,
                                      statistics_fields=None, multi_part="MULTI_PART", unsplit_lines="DISSOLVE_LINES"); _msg("10yr dissolve complete.")'''
        _msg(f"10yr Feature Class: {out_010}")

        _msg("All merges and dissolves complete.")

        # Cleanup temporaries
        _msg("Cleaning up temporary datasets...")
        for tmp in (merged_500, merged_100, merged_010):
            if arcpy.Exists(tmp): arcpy.management.Delete(tmp)
        _msg("Temporary datasets deleted.")

        # Derived outputs
        p_outfds.value = out_fds
        p_out500.value = out_500 if arcpy.Exists(out_500) else ""
        p_out100.value = out_100 if arcpy.Exists(out_100) else ""
        p_out10.value  = out_010 if arcpy.Exists(out_010) else ""

        _msg("Merge & Dissolve complete.")
        """
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
        p_ws  = arcpy.Parameter(
            displayName="Watershed",
            name="watershed",
            datatype="GPString",
            parameterType="Required",
            direction="Input"
        ); p_ws.filter.type = "ValueList"
        p_h8  = arcpy.Parameter(
            displayName="HUC8",
            name="huc8",
            datatype="GPString",
            parameterType="Required",
            direction="Input"
        ); p_h8.filter.type = "ValueList"
        p_gdb   = arcpy.Parameter(
            displayName="Project Geodatabase",
            name="project_gdb",
            datatype="Workspace",           # ← was DEWorkspace, Workspace is the safest for GDBs in a .pyt.
            parameterType="Required",
            direction="Input"
        )
        p_base  = arcpy.Parameter(
            displayName="Basename (from previous steps)",
            name="basename",
            datatype="GPString",
            parameterType="Required",
            direction="Input"
        )
        p_outfds= arcpy.Parameter(  # this should be updated to auto-populate based on previous selections made in the dialogue
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
        '''
        p_wsenv    = arcpy.Parameter(       # not used, envs applied globally
            displayName="Workspace",
            name="env_workspace",
            datatype="Workspace",           # ← was DEWorkspace, Workspace is the safest for GDBs in a .pyt.
            parameterType="Optional",
            direction="Input",
            category="Environment Settings"
        )
        p_scratch  = arcpy.Parameter(       # not used, envs applied globally
            displayName="Scratch Workspace",
            name="env_scratch",
            datatype="Workspace",           # ← was DEWorkspace, Workspace is the safest for GDBs in a .pyt.
            parameterType="Optional",
            direction="Input",
            category="Environment Settings"
        )
        p_ocs      = arcpy.Parameter(       # not used, envs applied globally
            displayName="Output Coordinate System",
            name="env_ocs",
            datatype="GPSpatialReference",
            parameterType="Optional",
            direction="Input",
            category="Environment Settings"
        )
'''
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

        return [p_lut, p_ws, p_h8, p_gdb, p_base, p_outfds, p_500, p_100, p_10,
                p_parallel, p_overwrite, #p_wsenv, p_scratch, p_ocs, # not used, envs applied globally
                p_validate,
                p_finalfds, p_topo, p_merged]

    def updateParameters(self, params):
        (p_lut, p_ws, p_h8, *_rest) = params[:3+10]  # just to keep readable
        p_basename = params[4] # derived

        # Cascade
        if p_lut.altered and p_lut.value:
            lookup = build_ws_huc_lookup(p_lut.valueAsText)
            p_ws.filter.list = sorted(lookup.keys())
            if p_ws.value and p_ws.value in lookup:
                p_h8.filter.list = lookup[p_ws.value]
            else:
                p_h8.filter.list = []

        if p_ws.value and p_h8.value:
            p_basename.value = compute_basename(p_ws.value, p_h8.value)
            params[5].value = os.path.join(params[3].valueAsText, f"Output_{p_basename.value}")
            params[6].value = os.path.join(params[3].valueAsText, f"Existing_Future_{p_basename.value}_500yr")
            params[7].value = os.path.join(params[3].valueAsText, f"Existing_Future_{p_basename.value}_100yr")
            params[8].value = os.path.join(params[3].valueAsText, f"Existing_Future_{p_basename.value}_10yr")
        return

    def execute(self, params, messages):
        (p_lut, p_ws, p_h8, p_gdb, p_base, p_outfds, p_500, p_100, p_10,
         p_parallel, p_overwrite, #p_wsenv, p_scratch, p_ocs, # not used, envs applied globally
         p_validate, # NEW
         p_finalfds, p_topo, p_merged) = params
        

        apply_envs(p_parallel.valueAsText or "100%",
                   p_overwrite.value,
                   #p_wsenv.valueAsText, # not used, envs applied globally
                   #p_scratch.valueAsText, # not used, envs applied globally
                   #p_ocs.value # not used, envs applied globally
                   )

        gdb       = p_gdb.valueAsText
        basename  = p_base.valueAsText
        out_fds   = p_outfds.valueAsText  # Output_<Basename> from Tool 2
        fc500     = p_500.valueAsText
        fc100     = p_100.valueAsText
        fc10      = p_10.valueAsText

        result = combine_floods(
            gdb=gdb,
            basename=basename,
            out_feature_dataset=out_fds,
            fc_500=fc500,
            fc_100=fc100,
            fc_10=fc10,
            #out_cs=p_ocs.value, # not used, envs applied globally
            validate_topology=bool(p_validate.value) # pass Flag
        )

        p_finalfds.value = result["final_fds"]
        p_topo.value     = result["topology"]
        p_merged.value   = result["ef_merged"]
        
        _msg("Combine (finalize & topology) completed.")

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
        p_ws_field = arcpy.Parameter( 
            displayName="Watershed Field",
            name="ws_field",
            datatype="Field",
            parameterType="Required",
            direction="Input"
        ); p_ws_field.parameterDependencies = [p_lut.name]
        
        p_huc_field = arcpy.Parameter(
            displayName="HUC Field",
            name="huc_field",
            datatype="Field",
            parameterType="Required",
            direction="Input"
        ); p_huc_field.parameterDependencies = [p_lut.name]
        '''
        p_ws  = arcpy.Parameter(
            displayName="Watershed",
            name="watershed",
            datatype="GPString",
            parameterType="Required",
            direction="Input"
        ); p_ws.filter.type="ValueList"
        p_h8  = arcpy.Parameter(
            displayName="HUC8",
            name="huc8",
            datatype="GPString",
            parameterType="Required",
            direction="Input"
        ); p_h8.filter.type="ValueList"

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
        '''
        p_table = arcpy.Parameter(      # not used, paths derived internally
            displayName="Input Table (with FC paths for Merge)",
            name="merge_input_table",
            datatype="Table View",
            parameterType="Required",
            direction="Input"
        )
        p_field = arcpy.Parameter(      # not used, paths derived internally
            displayName="Path Field Name",
            name="merge_path_field",
            datatype="GPString",
            parameterType="Optional",
            direction="Input"
        ); p_field.value="FCPath"
        '''

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
        '''
        p_wsenv    = arcpy.Parameter(      # not used, envs applied globally
            displayName="Workspace",
            name="env_workspace",
            datatype="Workspace",
            parameterType="Optional",
            direction="Input",
            category="Environment Settings"
        )
        p_scratch  = arcpy.Parameter(       # not used, envs applied globally
            displayName="Scratch Workspace",
            name="env_scratch",
            datatype="Workspace",
            parameterType="Optional",
            direction="Input",
            category="Environment Settings"
        )
        p_ocs      = arcpy.Parameter(       # not used, envs applied globally
            displayName="Output Coordinate System",
            name="env_ocs",
            datatype="GPSpatialReference",
            parameterType="Optional",
            direction="Input",
            category="Environment Settings"
        )
        '''

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

        return [p_lut, 
                #p_ws_field, p_huc_field, # Future enhancement
                p_ws, p_h8, p_gdb, p_sr,
                p_in_exist, p_in_future, p_splitfield,
                #p_table, p_field, # not used, paths derived internally
                p_parallel, p_overwrite, 
                #p_wsenv, p_scratch, p_ocs, # not used, envs applied globally
                p_basename, p_finalfds, p_topology, p_efmerged,
                p_validate_topo, # NEW
                p_skip_split, p_skip_merge, p_skip_combine]

    def updateParameters(self, params):
        (p_lut, p_ws, p_h8, *_rest) = params[:3+7]  # just to keep readable
        p_basename = params[10] # derived

        # Cascade
        if p_lut.altered and p_lut.value:
            lookup = build_ws_huc_lookup(p_lut.valueAsText)
            p_ws.filter.list = sorted(lookup.keys())
            if p_ws.value and p_ws.value in lookup:
                p_h8.filter.list = lookup[p_ws.value]
            else:
                p_h8.filter.list = []

        if p_ws.value and p_h8.value:
            p_basename.value = compute_basename(p_ws.value, p_h8.value)
        return
    
    ''' # Future enhancement: allow user to specify watershed and HUC fields if needed
    def updateParameters(self, params):
        # Extract parameters for readability
        p_lut, p_ws, p_h8, *_rest = params[:3+16]
        p_basename = params[12]  # Assuming this is a derived parameter

        if p_lut.altered and p_lut.value:
            # Build lookup to populate dependent lists
            lookup = build_ws_huc_lookup(p_lut.valueAsText, p_ws.valueAsText, p_h8.valueAsText)
            p_ws.filter.list = sorted(lookup.keys())

            if p_ws.value and p_ws.value in lookup:
                p_h8.filter.list = lookup[p_ws.value]
            else:
                p_h8.filter.list = []

        if p_ws.value and p_h8.value:
            p_basename.value = compute_basename(p_ws.value, p_h8.value)

        return
    '''

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
        (p_lut, 
         #p_ws_field, p_huc_field, 
         p_ws, p_h8, p_gdb, p_sr,
         p_in_exist, p_in_future, p_splitfield,
         #p_table, p_field,
         p_parallel, p_overwrite, 
         #p_wsenv, p_scratch, p_ocs,
         p_basename, p_finalfds, p_topology, p_efmerged,
         p_validate_topo, # NEW
         p_skip_split, p_skip_merge, p_skip_combine) = params

        # Apply envs once for all steps
        apply_envs(p_parallel.valueAsText or "100%",
                   p_overwrite.value,
                   #p_wsenv.valueAsText,
                   #p_scratch.valueAsText,
                   #p_ocs.value
                   )

        # Gather basics
        gdb         = p_gdb.valueAsText
        sr_obj      = p_sr.value
        existing_fc = p_in_exist.valueAsText
        future_fc   = p_in_future.valueAsText
        split_field = (p_splitfield.valueAsText or "FLOOD_FREQ")
        '''
        table_for_merge = p_table.valueAsText or os.path.join(gdb, "FeatureClassPaths") # Explicitly passed
        table_path_field= p_field.valueAsText or "FCPath" # Explicitly passed
        '''

        # Derived basename
        basename = compute_basename(p_ws.valueAsText, p_h8.valueAsText)
        arcpy.SetParameterAsText(12, basename)

        # Validate key inputs
        if not arcpy.Exists(gdb):
            raise RuntimeError(f"Project GDB not found: {gdb}")
        if not arcpy.Exists(existing_fc) or not arcpy.Exists(future_fc):
            raise RuntimeError("Existing/Future flood feature classes are required and must exist.")
        self._validate_field(existing_fc, split_field)
        self._validate_field(future_fc, split_field)

        # ---- 1) SPLIT ----
        if not p_skip_split.value:
            arcpy.AddMessage("1) SPLIT MODULE: Creating existing and future datasets and splitting by frequency…")
            _msg(f"Inputs:\n  Geodatabase: {gdb}\n  Existing FC: {existing_fc}\n  Future FC: {future_fc}\n  Split Field: {split_field}")
            split_out = split_floods(
                gdb=gdb,
                basename=basename,
                spatial_ref=sr_obj,
                existing_fc=existing_fc,
                future_fc=future_fc,
                split_field=split_field
            )
        else:
            arcpy.AddMessage("1) SPLIT MODULE: skipped by user flag.")
            split_out = {
                "fds_exist":  os.path.join(gdb, f"Existing_{basename}"),
                "fds_future": os.path.join(gdb, f"Future_{basename}"),
                "exist_10":   os.path.join(gdb, f"Existing_{basename}", f"Existing_{basename}_10yr"),
                "exist_100":  os.path.join(gdb, f"Existing_{basename}", f"Existing_{basename}_100yr"),
                "exist_500":  os.path.join(gdb, f"Existing_{basename}", f"Existing_{basename}_500yr"),
                "future_10":  os.path.join(gdb, f"Future_{basename}",  f"Future_{basename}_10yr"),
                "future_100": os.path.join(gdb, f"Future_{basename}",  f"Future_{basename}_100yr"),
                "future_500": os.path.join(gdb, f"Future_{basename}",  f"Future_{basename}_500yr"),
            }

        # Paths to feature classes that were modified
        fc_paths = {
            "Existing_10yr": split_out["exist_10"],
            "Existing_100yr": split_out["exist_100"],
            "Existing_500yr": split_out["exist_500"],
            "Future_10yr": split_out["future_10"],
            "Future_100yr": split_out["future_100"],
            "Future_500yr": split_out["future_500"],
        }
        _msg("Split step complete. Recording feature class paths to table for Merge step.")
        _msg(f"Feature classes created/modified:\n" + "\n".join(f"  {k}: {v}" for k,v in fc_paths.items()))
        
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

        arcpy.AddMessage(f"Feature class paths have been recorded in {output_table}.")
        #p_table.value = output_table  # set for potential reuse

        # ---- 2) MERGE ----
        # Reuse helpers from MergeFloodFreqs to keep behavior identical
        if not p_skip_merge.value:
            arcpy.AddMessage("2) MERGE SUBMODULE: Merging and dissolving by FLOOD_FREQ into output dataset.")
            merger = MergeFloodFreqs()
            paths = merger.read_fc_paths_from_table(output_table, "FCPath")
            if len(paths) < 3: raise RuntimeError("MERGE ERROR: expected at least 3 rows with FC paths in the input table.")
            buckets = merger.group_by_suffix(paths)
            sr_merge = merger.pick_spatial_ref(paths, sr_obj)
            out_fds  = merger.create_feature_dataset(gdb, basename, sr_merge)

            ''' # Not implemented, in_memory for merge temps where possible to reduce I/O
            _msg("Merging feature classes by frequency group into in_memory workspace...")
            tmp_500, tmp_100, tmp_010 = "in_memory\\m500", "in_memory\\m100", "in_memory\\m010"
            if buckets.get("500yr"): arcpy.management.Merge(buckets["500yr"], tmp_500)
            if buckets.get("100yr"): arcpy.management.Merge(buckets["100yr"], tmp_100)
            if buckets.get("10yr"):  arcpy.management.Merge(buckets["10yr"],  tmp_010)
            _msg("Merges complete.")
            '''
            
            # Temporary merged (on disk, for stability)
            merged_500 = os.path.join(out_fds, f"Merged_{basename}_500yr")
            merged_100 = os.path.join(out_fds, f"Merged_{basename}_100yr")
            merged_010 = os.path.join(out_fds, f"Merged_{basename}_10yr")
            
            # Final
            out_500 = os.path.join(out_fds, f"Existing_Future_{basename}_500yr")
            out_100 = os.path.join(out_fds, f"Existing_Future_{basename}_100yr")
            out_010 = os.path.join(out_fds, f"Existing_Future_{basename}_10yr")

            ''' # Not implemented, using disk-based temporaries for initial stability
            _msg("Dissolving merged feature classes...")
            if arcpy.Exists(tmp_500): merger.dissolve_on_freq(tmp_500, out_500, merger.DISSOLVE_FIELD)
            _msg(f"500yr dissolve complete: \n {out_500}")
            if arcpy.Exists(tmp_100): merger.dissolve_on_freq(tmp_100, out_100, merger.DISSOLVE_FIELD)
            _msg(f"100yr dissolve complete: \n {out_100}")
            if arcpy.Exists(tmp_010): merger.dissolve_on_freq(tmp_010, out_010, merger.DISSOLVE_FIELD)
            _msg(f"10yr dissolve complete: \n {out_010}")
            '''

            # Merge and Dissolve (simple, without normalization or safe fallback)
            _msg("Merging and dissolving feature classes by frequency group...")
            _msg("Merging 500yr group...")
            if buckets.get("500yr"):
                arcpy.management.Merge(buckets["500yr"] + buckets.get("100yr", []) + buckets.get("10yr", []), merged_500); _msg("500yr merge complete.")
            _msg("Repairing geometry on merged 500yr...")
            if arcpy.Exists(merged_500):
                try:
                    arcpy.management.RepairGeometry(merged_500, validation_method="OGC", delete_null="DELETE_NULL"); _msg("500yr repair complete.")
                except TypeError:
                    arcpy.management.RepairGeometry(merged_500, validation_method="OGC"); _msg("500yr repair complete (without delete_null).")
            if arcpy.Exists(merged_500):
                try:
                    # Add a new field for flood frequency if it doesn't exist
                    fields = [field.name for field in arcpy.ListFields(merged_500)]
                    if merger.DISSOLVE_FIELD not in fields:
                        arcpy.management.AddField(merged_500, merger.DISSOLVE_FIELD, "Text")
                    
                    # Calculate the flood frequency value as 0.2
                    arcpy.management.CalculateField(merged_500, merger.DISSOLVE_FIELD, "0.2", "PYTHON3")
                    _msg("Flood frequency calculation complete.")
                except Exception as e:
                    _msg(f"Error during flood frequency calculation: {e}")
            _msg("Dissolving 500yr group...")
            if arcpy.Exists(merged_500):
                arcpy.management.Dissolve(merged_500, out_500, dissolve_field=merger.DISSOLVE_FIELD,
                                        statistics_fields=None, multi_part="MULTI_PART", unsplit_lines="DISSOLVE_LINES"); _msg("500yr dissolve complete.")
            _msg(f"500yr Feature Class: {out_500}")
            
            _msg("Merging 100yr group...")
            if buckets.get("100yr"):
                arcpy.management.Merge(buckets["100yr"] + buckets.get("10yr", []), merged_100); _msg("100yr merge complete.")
            _msg("Repairing geometry on merged 100yr...")
            if arcpy.Exists(merged_100):
                try:
                    arcpy.management.RepairGeometry(merged_100, validation_method="OGC", delete_null="DELETE_NULL"); _msg("100yr repair complete.")
                except TypeError:
                    arcpy.management.RepairGeometry(merged_100, validation_method="OGC"); _msg("100yr repair complete (without delete_null).")
            if arcpy.Exists(merged_100):
                try:
                    # Add a new field for flood frequency if it doesn't exist
                    fields = [field.name for field in arcpy.ListFields(merged_100)]
                    if merger.DISSOLVE_FIELD not in fields:
                        arcpy.management.AddField(merged_100, merger.DISSOLVE_FIELD, "Text")
                    
                    # Calculate the flood frequency value as 1
                    arcpy.management.CalculateField(merged_100, merger.DISSOLVE_FIELD, "1", "PYTHON3")
                    _msg("Flood frequency calculation complete.")
                except Exception as e:
                    _msg(f"Error during flood frequency calculation: {e}")
            _msg("Dissolving 100yr group...")
            if arcpy.Exists(merged_100):
                arcpy.management.Dissolve(merged_100, out_100, dissolve_field=merger.DISSOLVE_FIELD,
                                        statistics_fields=None, multi_part="MULTI_PART", unsplit_lines="DISSOLVE_LINES"); _msg("100yr dissolve complete.")
            _msg(f"100yr Feature Class: {out_100}")
            
            _msg("Merging 10yr group...")
            if buckets.get("10yr"):
                arcpy.management.Merge(buckets["10yr"], merged_010); _msg("10yr merge complete.")
            _msg("Repairing geometry on merged 10yr...")
            if arcpy.Exists(merged_010):
                try:
                    arcpy.management.RepairGeometry(merged_010, validation_method="OGC", delete_null="DELETE_NULL"); _msg("10yr repair complete.")
                except TypeError:
                    arcpy.management.RepairGeometry(merged_010, validation_method="OGC"); _msg("10yr repair complete (without delete_null).")
            if arcpy.Exists(merged_010):
                try:
                    # Add a new field for flood frequency if it doesn't exist
                    fields = [field.name for field in arcpy.ListFields(merged_010)]
                    if merger.DISSOLVE_FIELD not in fields:
                        arcpy.management.AddField(merged_010, merger.DISSOLVE_FIELD, "Text")
                    
                    # Calculate the flood frequency value as 10
                    arcpy.management.CalculateField(merged_010, merger.DISSOLVE_FIELD, "10", "PYTHON3")
                    _msg("Flood frequency calculation complete.")
                except Exception as e:
                    _msg(f"Error during flood frequency calculation: {e}")
            _msg("Dissolving 10yr group...")
            if arcpy.Exists(merged_010):
                arcpy.management.Dissolve(merged_010, out_010, dissolve_field=merger.DISSOLVE_FIELD,
                                        statistics_fields=None, multi_part="MULTI_PART", unsplit_lines="DISSOLVE_LINES"); _msg("10yr dissolve complete.")
            _msg(f"10yr Feature Class: {out_010}")

            '''# Cleanup local temp files
            for tmp in (merged_500, merged_100, merged_010): # Leaving temps for now for debugging
                if arcpy.Exists(tmp):
                    arcpy.management.Delete(tmp)'''

            '''
            # Cleanup in_memory temp files
            _msg("Cleaning up temporary datasets...")
            for t in (tmp_500, tmp_100, tmp_010):
                if arcpy.Exists(t): arcpy.management.Delete(t)
            _msg("Temporary datasets deleted.")
            '''
        else:
            arcpy.AddMessage("2) MERGE SUBMODULE: skipped by user flag.")
            out_fds = os.path.join(gdb, f"Output_{basename}")
            out_500 = os.path.join(out_fds, f"Existing_Future_{basename}_500yr")
            out_100 = os.path.join(out_fds, f"Existing_Future_{basename}_100yr")
            out_010 = os.path.join(out_fds, f"Existing_Future_{basename}_10yr")
            _msg(f"Assuming existing output dataset and feature classes at: \n{out_500}, \n{out_100}, \n{out_010}")
            if not any(arcpy.Exists(fc) for fc in (out_500, out_100, out_010)):
                raise RuntimeError("WARNING: MERGE SUBMODULE skipped but no existing flood layers found in output dataset.")

        # ---- 3) COMBINE ----
        if not p_skip_combine.value:
            arcpy.AddMessage("3) COMBINE SUBMODULE: Repair, erase overlaps, merge, and (optionally) validate topology…")
            _msg(f"Inputs:\n  Flood_500yr: {out_500}\n  Flood_100yr: {out_100}\n  Flood_10yr: {out_010}")
            combine_out = combine_floods(
                gdb=gdb,
                basename=basename,
                out_feature_dataset=out_fds,   # from Merge section
                fc_500=out_500,
                fc_100=out_100,
                fc_10=out_010,
                #out_cs=p_ocs.value,
                validate_topology=bool(p_validate_topo.value)  # ← pass checkbox
            )
            p_finalfds.value = combine_out["final_fds"]
            p_topology.value = combine_out["topology"]
            p_efmerged.value = combine_out["ef_merged"]
        else:
            arcpy.AddMessage("3) COMBINE SUBMODULE: skipped by user flag.")

        arcpy.AddMessage("Flood Quilt pipeline complete.")
