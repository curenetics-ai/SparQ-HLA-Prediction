#!/bin/bash
# Entrypoint for the NetChop container. Needs bash for process substitution.
#
#   selftest          verify the installation (see checks below)
#   <anything else>   executed verbatim, e.g. `netChop -t 0.5 -v 1 in.fsa`
#
# NOTE ON THE VENDOR FIXTURE
# The distribution ships test/test.fsa and test/test.out, but they do not
# correspond: test.fsa contains Influenza A M2 (BHB191648, MSLLTEVETPI...)
# while test.out holds predictions for gi|3333147 (M, A, G...), generated in
# 2006 from a different input. Diffing them therefore proves nothing, so the
# selftest below checks properties that must hold for any correct install.
set -eu

if [ "${1:-selftest}" = "selftest" ]; then
    fsa="$NETCHOP/test/test.fsa"
    fail=0
    note() { printf '  %-6s %s\n' "$1" "$2"; }

    echo "NETCHOP=$NETCHOP"
    echo "networks available: $(ls "$NETCHOP/etc" | tr '\n' ' ')"
    echo

    cterm=$(mktemp); s20=$(mktemp); repeat=$(mktemp)
    netChop -t 0.500000 -v 0 "$fsa" > "$cterm"
    netChop -t 0.500000 -v 1 "$fsa" > "$s20"
    netChop -t 0.500000 -v 0 "$fsa" > "$repeat"

    rows() { grep -E '^ *[0-9]' "$1"; }
    scores() { rows "$1" | awk '{print $4}'; }

    # 1. the binary produced predictions at all
    n=$(rows "$cterm" | wc -l)
    if [ "$n" -gt 0 ]; then note OK "produced $n prediction rows"
    else note FAIL "no prediction rows"; fail=1; fi

    # 2. the two networks are genuinely different models, not an ignored flag.
    #    This is the check that matters: selecting 20S must actually change the
    #    output, otherwise -v is silently a no-op and every result is C-term.
    if [ "$(scores "$cterm")" = "$(scores "$s20")" ]; then
        note FAIL "-v 0 and -v 1 gave identical scores: network flag is a no-op"
        fail=1
    else
        d=$(paste -d' ' <(scores "$cterm") <(scores "$s20") \
            | awk '$1!=$2{c++} END{print c+0}')
        note OK "-v 0 (Cterm) and -v 1 (20S) differ at $d/$n positions"
    fi

    # 3. the header must name the network actually used
    if grep -q 'version C-term' "$cterm"; then note OK "-v 0 header reports C-term"
    else note FAIL "-v 0 header does not report C-term"; fail=1; fi
    if grep -qi 'version 20S' "$s20"; then note OK "-v 1 header reports 20S"
    else note FAIL "-v 1 header does not report 20S"; fail=1; fi

    # 4. determinism - same input, same output
    if [ "$(scores "$cterm")" = "$(scores "$repeat")" ]; then
        note OK "deterministic across runs"
    else note FAIL "non-deterministic output"; fail=1; fi

    # 5. scores are probabilities, and the S/. flag agrees with the threshold
    bad=$(scores "$s20" | awk '$1<0 || $1>1 {c++} END{print c+0}')
    if [ "$bad" -eq 0 ]; then note OK "all scores within [0,1]"
    else note FAIL "$bad scores outside [0,1]"; fail=1; fi

    mism=$(rows "$s20" | awk '{s=($3=="S"); c=($4>=0.5); if (s!=c) n++} END{print n+0}')
    if [ "$mism" -eq 0 ]; then note OK "cleavage flag consistent with threshold 0.5"
    else note FAIL "$mism rows where S/. disagrees with the threshold"; fail=1; fi

    echo
    if [ "$fail" -eq 0 ]; then echo "SELFTEST PASS"; exit 0; fi
    echo "SELFTEST FAIL"; exit 1
fi

exec "$@"
