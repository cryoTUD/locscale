# vim: set expandtab shiftwidth=4 softtabstop=4:
"""
Delft University of Technology (TU Delft) hereby disclaims all copyright interest in the
program 'LocScale2' written by the Author(s).

Copyright (C) 2026 Alok Bharadwaj and Arjen J. Jakobi
"""

from chimerax.core.toolshed import BundleAPI


class _LocScale2API(BundleAPI):

    api_version = 1

    @staticmethod
    def start_tool(session, bi, ti):
        from . import tool
        return tool.LocScale2Tool(session, ti.name)

    @staticmethod
    def register_command(bi, ci, logger):
        from chimerax.core.commands import register
        if ci.name == "locscale2 verify":
            from . import verify
            register(ci.name, verify.verify_desc, verify.verify, logger=logger)
        else:
            from . import cmd
            register(ci.name, cmd.locscale2_desc, cmd.locscale2, logger=logger)


bundle_api = _LocScale2API()
