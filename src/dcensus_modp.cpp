// Level-by-level computation of the down-census vectors d_j(lambda) modulo a
// 61/62-bit prime for ALL partitions of n, n = 0..N, with collision detection.
//
// If d(lambda) = d(mu) exactly then d(lambda) == d(mu) (mod P); hence "no
// collision mod P" rigorously implies "no exact collision".  Any collision
// found mod P must be re-checked exactly (Python).
//
// Output per level: n p(n) #distinct #symmetric #transpose_pairs #collision_groups digest
// where digest = sum_{lambda, j} d_j(lambda) * 7^j  (mod P)  [order independent;
// the same digest is printed by the exact Python scan for cross-checking].
//
// Build: g++ -O2 -std=c++17 -o dcensus_modp dcensus_modp.cpp
// Usage: ./dcensus_modp N [prime_index]   (prime_index 0: 2^61-1, 1: 2^62-57)
#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <vector>

typedef unsigned long long u64;
typedef unsigned __int128 u128;

static u64 P = (1ULL << 61) - 1;
static const int MAXN = 90;
static u64 Ptab[MAXN + 2][MAXN + 2];  // Ptab[m][k] = #partitions of m with parts <= k

static int N;
static int curN;
static std::vector<u64> prevv, curv;       // rows of length curN (prev) / curN+1 (cur)
static std::vector<uint8_t> parts_flat;    // partitions of level curN, concatenated
static std::vector<uint32_t> parts_off;    // offsets (size p(n)+1)
static u64 counter;
static uint8_t parts[MAXN + 2];

static inline u64 addmod(u64 a, u64 b) { u64 c = a + b; return c >= P ? c - P : c; }

static u64 rank_partition(const uint8_t* q, int len, int n) {
    u64 r = 0; int rem = n; int prev = n;
    for (int i = 0; i < len; ++i) {
        int hi = std::min(prev, rem);
        for (int s = q[i] + 1; s <= hi; ++s) r += Ptab[rem - s][s];
        rem -= q[i]; prev = q[i];
    }
    return r;
}

static void process(int len) {
    // parts[0..len) is a partition of curN, in reverse-lex order, index = counter
    u64 idx = counter++;
    int n = curN;
    u64* row = &curv[idx * (n + 1)];
    row[0] = 1;
    for (int j = 1; j <= n; ++j) row[j] = 0;
    uint8_t q[MAXN + 2];
    for (int i = 0; i < len; ++i) q[i] = parts[i];
    for (int i = 0; i < len; ++i) {
        int nxt = (i + 1 < len) ? parts[i + 1] : 0;
        if (parts[i] > nxt) {  // removable corner in row i
            q[i] = parts[i] - 1;
            int qlen = (q[i] == 0) ? len - 1 : len;
            u64 ridx = rank_partition(q, qlen, n - 1);
            const u64* prow = &prevv[ridx * n];
            for (int j = 1; j <= n; ++j) row[j] = addmod(row[j], prow[j - 1]);
            q[i] = parts[i];
        }
    }
    parts_off.push_back((uint32_t)parts_flat.size());
    for (int i = 0; i < len; ++i) parts_flat.push_back(parts[i]);
}

static void gen(int rem, int maxpart, int len) {
    if (rem == 0) { process(len); return; }
    for (int s = std::min(rem, maxpart); s >= 1; --s) {
        parts[len] = (uint8_t)s;
        gen(rem - s, s, len + 1);
    }
}

static inline u64 mix(u64 x) {
    x += 0x9E3779B97F4A7C15ULL;
    x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
    x = (x ^ (x >> 27)) * 0x94D049BB133111EBULL;
    return x ^ (x >> 31);
}

static void conjugate(const uint8_t* q, int len, std::vector<uint8_t>& out) {
    out.clear();
    if (len == 0) return;
    for (int j = 0; j < q[0]; ++j) {
        int c = 0;
        for (int i = 0; i < len; ++i) if (q[i] > j) ++c;
        out.push_back((uint8_t)c);
    }
}

static void print_partition(FILE* f, const uint8_t* q, int len) {
    fprintf(f, "(");
    for (int i = 0; i < len; ++i) fprintf(f, "%d%s", (int)q[i], i + 1 < len ? "," : "");
    fprintf(f, ")");
}

int main(int argc, char** argv) {
    N = argc > 1 ? atoi(argv[1]) : 60;
    int pidx = argc > 2 ? atoi(argv[2]) : 0;
    if (pidx == 1) P = 4611686018427387847ULL;  // 2^62 - 57, prime
    if (N > MAXN) { fprintf(stderr, "N too large\n"); return 1; }
    // partition-count table
    for (int m = 0; m <= MAXN + 1; ++m)
        for (int k = 0; k <= MAXN + 1; ++k) {
            if (m == 0) Ptab[m][k] = 1;
            else if (k == 0) Ptab[m][k] = 0;
            else if (k > m) Ptab[m][k] = Ptab[m][m];
            else Ptab[m][k] = Ptab[m][k - 1] + Ptab[m - k][k];
        }
    printf("# P = %llu\n# n p(n) distinct symmetric pairs collisions digest\n", P);
    prevv.assign(1, 1);  // level 0: the empty partition, d = (1)
    {
        // report level 0
        printf("0 1 1 1 0 0 1\n"); fflush(stdout);
    }
    std::vector<uint8_t> conj;
    for (int n = 1; n <= N; ++n) {
        curN = n;
        u64 pn = Ptab[n][n];
        curv.assign(pn * (n + 1), 0);
        parts_flat.clear(); parts_off.clear();
        counter = 0;
        gen(n, n, 0);
        if (counter != pn) { fprintf(stderr, "count mismatch at n=%d\n", n); return 1; }
        parts_off.push_back((uint32_t)parts_flat.size());
        // digest
        u64 dig = 0; u64 w = 1;
        std::vector<u64> ws(n + 1);
        for (int j = 0; j <= n; ++j) { ws[j] = w; w = (u64)(((u128)w * 7) % P); }
        for (u64 i = 0; i < pn; ++i) {
            const u64* row = &curv[i * (n + 1)];
            u128 acc = 0;
            for (int j = 0; j <= n; ++j) acc += (u128)row[j] * ws[j];
            dig = addmod(dig, (u64)(acc % P));
        }
        // hashes and collision detection
        std::vector<std::pair<u64, u64>> hv(pn);
        for (u64 i = 0; i < pn; ++i) {
            const u64* row = &curv[i * (n + 1)];
            u64 h = 0x1234567ULL;
            for (int j = 0; j <= n; ++j) h = mix(h ^ row[j]);
            hv[i] = {h, i};
        }
        std::sort(hv.begin(), hv.end());
        u64 distinct = 0, nsym = 0, npairs = 0, ncoll = 0;
        u64 a = 0;
        while (a < pn) {
            u64 b = a;
            while (b < pn && hv[b].first == hv[a].first) ++b;
            // within [a,b): group by exact row equality
            std::vector<char> used(b - a, 0);
            for (u64 x = a; x < b; ++x) {
                if (used[x - a]) continue;
                std::vector<u64> grp; grp.push_back(hv[x].second); used[x - a] = 1;
                const u64* rx = &curv[hv[x].second * (n + 1)];
                for (u64 y = x + 1; y < b; ++y) {
                    if (used[y - a]) continue;
                    const u64* ry = &curv[hv[y].second * (n + 1)];
                    bool eq = true;
                    for (int j = 0; j <= n; ++j) if (rx[j] != ry[j]) { eq = false; break; }
                    if (eq) { grp.push_back(hv[y].second); used[y - a] = 1; }
                }
                ++distinct;
                // classify group
                const uint8_t* q0 = &parts_flat[parts_off[grp[0]]];
                int l0 = parts_off[grp[0] + 1] - parts_off[grp[0]];
                conjugate(q0, l0, conj);
                bool sym = ((int)conj.size() == l0) && std::equal(conj.begin(), conj.end(), q0);
                bool ok = false;
                if (grp.size() == 1 && sym) { ++nsym; ok = true; }
                else if (grp.size() == 2 && !sym) {
                    const uint8_t* q1 = &parts_flat[parts_off[grp[1]]];
                    int l1 = parts_off[grp[1] + 1] - parts_off[grp[1]];
                    if ((int)conj.size() == l1 && std::equal(conj.begin(), conj.end(), q1)) { ++npairs; ok = true; }
                }
                if (!ok) {
                    ++ncoll;
                    printf("COLLISION n=%d:", n);
                    for (u64 g : grp) { printf(" "); print_partition(stdout, &parts_flat[parts_off[g]], parts_off[g + 1] - parts_off[g]); }
                    printf("\n");
                }
            }
            a = b;
        }
        printf("%d %llu %llu %llu %llu %llu %llu\n", n, pn, distinct, nsym, npairs, ncoll, dig);
        fflush(stdout);
        prevv.swap(curv);
        std::vector<u64>().swap(curv);
    }
    return 0;
}
