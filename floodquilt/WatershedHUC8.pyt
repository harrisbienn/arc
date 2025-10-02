# filename: WatershedHUC8.pyt
import arcpy
import re
from collections import defaultdict

def sanitize(s):
    """Remove whitespace and special characters, keep only letters, numbers, and underscores."""
    return re.sub(r'[^A-Za-z0-9]', '', s)

class Toolbox(object):
    def __init__(self):
        self.label = "Watershed/HUC8 Toolbox"
        self.alias = "ws_huc8"
        self.tools = [PickWatershedHUC8]

class PickWatershedHUC8(object):
    def __init__(self):
        self.label = "Pick Watershed & HUC8"
        self.description = "Select Watershed and HUC8 from a table with multiple columns."
        self.canRunInBackground = False

    def getParameterInfo(self):
        table = arcpy.Parameter(
            displayName="Lookup Table (two columns: Watershed, HUC8)",
            name="table",
            datatype="DETable",
            parameterType="Required",
            direction="Input"
        )

        ws = arcpy.Parameter(
            displayName="TWDB Watershed",
            name="TWDBwatershed",
            datatype="GPString",
            parameterType="Required",
            direction="Input"
        )
        ws.filter.type = "ValueList"

        huc8 = arcpy.Parameter(
            displayName="HUC8 Name",
            name="huc8",
            datatype="GPString",
            parameterType="Required",
            direction="Input"
        )
        huc8.filter.type = "ValueList"

        # Optional output path showing inline use
        out = arcpy.Parameter(
            displayName="Output name (uses {Watershed}_{HUC8})",
            name="out_name",
            datatype="GPString",
            parameterType="Optional",
            direction="Input"
        )
        return [table, ws, huc8, out]

    def updateParameters(self, params):
        table, ws, huc8, out = params

        if table.altered and table.value:
            ws_to_huc = defaultdict(set)
            fields = [f.name for f in arcpy.ListFields(table.value)]

            def findcol(want):
                for f in fields:
                    if f.lower() == want.lower():
                        return f
                return want

            ws_col = findcol("TWDBWatershed")
            huc_col = findcol("Name")

            with arcpy.da.SearchCursor(table.value, [ws_col, huc_col]) as cur:
                for w, h in cur:
                    if w is not None and h is not None:
                        ws_to_huc[str(w)].add(str(h))

            params[1].filter.list = sorted(ws_to_huc.keys())

            if ws.value and ws.value in ws_to_huc:
                params[2].filter.list = sorted(ws_to_huc[ws.value])
            else:
                params[2].filter.list = []

        # Auto-compose sanitized output name
        if ws.value and huc8.value and not out.altered:
            clean_ws = sanitize(ws.value)
            clean_huc8 = sanitize(huc8.value)
            out.value = f"{clean_ws}_{clean_huc8}"

        return

    def execute(self, params, messages):
        ws = params[1].valueAsText
        huc8 = params[2].valueAsText
        out = params[3].valueAsText if params[3].valueAsText else f"{ws}_{huc8}"

        # Always sanitize in final output
        clean_ws = sanitize(ws)
        clean_huc8 = sanitize(huc8)
        clean_out = sanitize(out) if "_" not in out else f"{clean_ws}_{clean_huc8}"

        messages.addMessage(f"Selected Watershed={ws}, HUC8={huc8}, OutputName={clean_out}")
        return