__version__ = "2.0"

import contextlib
from typing import Optional, List, TextIO

from gecos.gecos_rdkit import GecosRdkit
from gecos.gecos_pybabel import GecosPyBabel
from gecos.gecos_extract_neighbors import GecosExtractNeighbors
from gecos.send_qm_conformers import send_qm_conformers, check_qm_jobs


def main(argv: Optional[List[str]] = None, stream: Optional[TextIO] = None) -> int:

    from gecos_analysis.gecos_analysis import main_app

    with contextlib.ExitStack() as ctx:
        return main_app()
