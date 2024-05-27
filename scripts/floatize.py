#!/usr/bin/env python3

import pathlib
import re

DOUBLE_LITERAL_REGEX = re.compile(r"\d+\.\d+(?![ef])", re.IGNORECASE)
EXP_DOUBLE_LITERAL_REGEX = re.compile(r"(\d+e[+-]?)(\d+)(?!f)", re.IGNORECASE)


rootDir = pathlib.Path(__file__).parent.parent  # root directory of the project

incDir = rootDir / "include" / "ruckig"  # include directory
srcDir = rootDir / "src" / "ruckig"  # source directory


for dir in [incDir, srcDir]:
    for file in dir.iterdir():
        if file.suffix in {".h", ".hpp", ".cpp"}:
            with file.open("r+") as f:
                content = f.read()

                content = content.replace("double", "float")
                content = content.replace("DBL_EPSILON", "FLT_EPSILON")
                content = DOUBLE_LITERAL_REGEX.sub(r"\g<0>f", content)
                content = EXP_DOUBLE_LITERAL_REGEX.sub(
                    lambda m: f"{m[1]}{int(m[2])//(2 if int(m[2])>7 else 1)}f", content
                )

                f.seek(0)
                f.write(content)
                f.truncate()
