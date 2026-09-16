// topend3_scan: for each n in [n1, n2], enumerate all partitions of n, compute the content
// moments C2, C1^2, C4 (exact int64) and log f^lam (double), canonicalise under transpose,
// and print every pair of distinct transpose classes with equal (C2, C1^2, C4) and
// |log f - log f'| < 1e-7  (candidates for equal (n, f, u_3, u_4, u_5); exact confirmation in Python).
// usage: topend3_scan n1 n2
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <algorithm>
#include <cstdint>
using namespace std;
typedef unsigned __int128 u128;
struct Rec { int64_t C4; int32_t C2; int32_t C1a; double logf; u128 enc; int len; };
static int n;
static vector<Rec> recs;
static vector<int> parts;
static double logfact[200];
static double lg[400];

static void process() {
    int l = parts.size();
    // conjugate
    int l1 = parts[0];
    static int conj[200];
    for (int j = 0; j < l1; j++) { int c = 0; for (int i = 0; i < l; i++) if (parts[i] > j) c++; conj[j] = c; }
    int64_t c1 = 0, c2 = 0, c4 = 0;
    double lf = logfact[n];
    for (int i = 0; i < l; i++) for (int j = 0; j < parts[i]; j++) {
        int64_t c = j - i; c1 += c; c2 += c * c; c4 += c * c * c * c;
        int h = parts[i] - j + conj[j] - i - 1; lf -= lg[h];
    }
    // boundary word: for i = l-1 .. 0: zeros (parts[i]-parts[i+1]) then a one -> read from bottom-left.
    // word W (length parts[0]+l): bits appended MSB first.
    u128 w = 0, wt = 0; int len = l1 + l;
    // build w
    {   int prev = 0;  // parts[l] = 0
        for (int i = l - 1; i >= 0; i--) { int z = parts[i] - prev; for (int t = 0; t < z; t++) w = (w << 1); w = (w << 1) | 1; prev = parts[i]; }
    }
    // transpose word = reverse + complement of w
    {   u128 x = w; for (int t = 0; t < len; t++) { wt = (wt << 1) | (((x >> t) & 1) ^ 1); } }
    u128 enc = (w < wt) ? w : wt;
    Rec r; r.len = len; r.C4 = c4; r.C2 = (int32_t)c2; r.C1a = (int32_t)(c1 < 0 ? -c1 : c1); r.logf = lf; r.enc = enc;
    recs.push_back(r);
}
static void gen(int rem, int maxp) {
    if (rem == 0) { process(); return; }
    for (int p = min(rem, maxp); p >= 1; p--) { parts.push_back(p); gen(rem - p, p); parts.pop_back(); }
}
static void decode(u128 enc, int len, vector<int>& out) {
    vector<int> p; int zeros = 0;
    for (int t = len - 1; t >= 0; t--) { int b = (int)((enc >> t) & 1); if (b == 0) zeros++; else p.push_back(zeros); }
    out.assign(p.rbegin(), p.rend());
}
int main(int argc, char** argv) {
    int n1 = atoi(argv[1]), n2 = atoi(argv[2]);
    logfact[0] = 0; for (int i = 1; i < 200; i++) logfact[i] = logfact[i - 1] + log((double)i);
    for (int i = 1; i < 400; i++) lg[i] = log((double)i);
    for (n = n1; n <= n2; n++) {
        recs.clear(); parts.clear(); gen(n, n);
        sort(recs.begin(), recs.end(), [](const Rec& a, const Rec& b) {
            if (a.C2 != b.C2) return a.C2 < b.C2; if (a.C1a != b.C1a) return a.C1a < b.C1a;
            if (a.C4 != b.C4) return a.C4 < b.C4; return a.logf < b.logf; });
        long groups5 = 0, cand = 0;
        size_t i = 0;
        while (i < recs.size()) {
            size_t j = i; while (j < recs.size() && recs[j].C2 == recs[i].C2 && recs[j].C1a == recs[i].C1a && recs[j].C4 == recs[i].C4) j++;
            // within [i,j): classes distinct enc; sorted by logf
            if (j - i >= 2) {
                vector<Rec> u;  // unique classes in the group, sorted by logf
                for (size_t a = i; a < j; a++) { bool dup = false; for (auto& q : u) if (q.enc == recs[a].enc) { dup = true; break; } if (!dup) u.push_back(recs[a]); }
                if (u.size() >= 2) groups5++;
                for (size_t a = 0; a < u.size(); a++) for (size_t b = a + 1; b < u.size() && u[b].logf - u[a].logf < 1e-7; b++) {
                    vector<int> pa, pb; decode(u[a].enc, u[a].len, pa); decode(u[b].enc, u[b].len, pb);
                    printf("CAND n=%d [", n); for (size_t t = 0; t < pa.size(); t++) printf("%d%s", pa[t], t + 1 < pa.size() ? "," : "");
                    printf("] ["); for (size_t t = 0; t < pb.size(); t++) printf("%d%s", pb[t], t + 1 < pb.size() ? "," : ""); printf("]\n");
                    cand++;
                }
            }
            i = j;
        }
        printf("SUMMARY n=%d partitions=%zu groups_with_equal_(C2,C1sq,C4)_and_>=2_classes=%ld candidate_pairs_equal_logf=%ld\n", n, recs.size(), groups5, cand);
        fflush(stdout);
    }
    return 0;
}
