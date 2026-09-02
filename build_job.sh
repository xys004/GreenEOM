#!/usr/bin/env bash
# Concatenate the package and a script into one self-contained job for the
# Mathematica bridge, redirecting Print into a log that is returned as the
# job's value (the bridge returns the last expression, not the stdout stream).
set -euo pipefail
here="$(cd "$(dirname "$0")" && pwd)"
script="${1:-tests.wl}"
out="${2:-job.wl}"

{
  # The bridge kernel is persistent, so wipe any earlier copy of the package
  # first: otherwise every symbol comes back as a ::shdw warning and, worse,
  # a stale definition could be the one that answers.
  echo 'Quiet[Remove["GreenEOM`*", "GreenEOM`Private`*", "GreenNEQ`*", "GreenNEQ`Private`*"]];'
  echo ''
  cat "$here/GreenEOM.wl"
  echo ''
  cat "$here/GreenNEQ.wl"
  echo ''
  echo '$log = {};'
  echo 'lg[args___] := AppendTo[$log, StringJoin @@ (If[StringQ[#], #, ToString[#]] & /@ {args})];'
  echo ''
  sed 's/Print\[/lg[/g' "$here/$script"
  echo ''
  echo 'StringRiffle[$log, "\n"]'
} > "$here/$out"

echo "wrote $here/$out ($(wc -l < "$here/$out") lines)"
