from datetime import datetime
from typing import Any


def add_warning(warnings: list[dict[str, Any]], level: str, source: str, message: str, context: Any = None) -> None:
    warnings.append(
        {
            "Timestamp": datetime.now().isoformat(timespec="seconds"),
            "Level": level,
            "Source": source,
            "Message": message,
            "Context": context,
        }
    )
