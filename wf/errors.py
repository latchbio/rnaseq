from dataclasses import dataclass
from enum import Enum
from typing import Dict, Optional
from dataclasses_json import dataclass_json


class VerifiedLogLevel(str, Enum):
    WARNING = "warning"
    ERROR = "error"


@dataclass_json
@dataclass
class VerifiedLog(Exception):
    level: VerifiedLogLevel
    name: str
    step: str
    traceback: str
    description: str
    file_previews: Optional[Dict[str, str]]


class LatchValidationError(Exception):
    pass
