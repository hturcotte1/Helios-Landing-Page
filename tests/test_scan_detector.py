"""'Test the tester': the collision detectors must fire on a synthetic collision."""
import os, sys, subprocess
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from dcensus_scan import level_dvectors, collision_report
from young import partitions, conjugate


def test_python_detector_fires():
    prev = {(): (1,)}
    for n in range(1, 9):
        prev = level_dvectors(n, prev)
    groups, bad, n_sym, digest, dig_p = collision_report(8, prev)
    assert bad == [] and n_sym == 2 and len(groups) == 12
    fake = dict(prev)
    fake[(6, 2)] = fake[(2, 2, 1, 1, 1, 1)] = fake[(5, 3)]   # synthetic collision (both members of a transpose pair)
    groups, bad, n_sym, digest, dig_p = collision_report(8, fake)
    assert len(bad) == 1 and set(bad[0][1]) == {(6, 2), (2, 2, 1, 1, 1, 1), (5, 3), (2, 2, 2, 1, 1)}
    # overwriting only one member of a transpose pair leaves a lone non-symmetric shape, which is
    # flagged too (the report raises, since d(lam) != d(lam^t) is impossible for genuine data)
    fake = dict(prev); fake[(6, 2)] = fake[(5, 3)]
    import pytest
    with pytest.raises(AssertionError):
        collision_report(8, fake)


def test_cpp_detector_fires():
    src = os.path.join(os.path.dirname(__file__), "..", "src", "dcensus_modp.cpp")
    code = open(src).read()
    # inject: at level 8, copy the row of partition index 1 ((7,1)) over index 2 ((6,2)) before detection
    marker = "        // digest\n"
    assert marker in code
    inj = ("        if (n == 8) { for (int j = 0; j <= n; ++j) curv[2 * (n + 1) + j] = curv[1 * (n + 1) + j]; }\n")
    tmp = "/tmp/dcensus_modp_inject.cpp"
    open(tmp, "w").write(code.replace(marker, marker + inj))
    exe = "/tmp/dcensus_modp_inject"
    subprocess.run(["g++", "-O2", "-std=c++17", "-o", exe, tmp], check=True)
    out = subprocess.run([exe, "9"], capture_output=True, text=True, check=True).stdout
    coll = [l for l in out.splitlines() if l.startswith("COLLISION n=8")]
    # two flagged groups: {(7,1),(2,1^6),(6,2)} and the orphaned singleton (2,2,1,1,1,1)
    assert len(coll) == 2 and any("(7,1)" in l and "(6,2)" in l for l in coll), out
    assert any("(2,2,1,1,1,1)" in l for l in coll), out
    lines = [l for l in out.splitlines() if l.startswith("8 ")]
    assert lines and lines[0].split()[5] == "2"
