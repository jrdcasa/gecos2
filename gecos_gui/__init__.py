import contextlib
from typing import Optional, List, TextIO


def main(argv: Optional[List[str]] = None, stream: Optional[TextIO] = None) -> int:

    from gecos_gui.maingui import main_gui_app

    with contextlib.ExitStack() as ctx:
        return main_gui_app()
