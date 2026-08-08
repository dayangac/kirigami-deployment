// freeze_corpus.cpp -- writes every input population the kill / yield / regime / scaling
// experiments and the tests build procedurally to JSON, so that a port (kirigami-julia)
// can rerun the experiments on byte-identical inputs without re-implementing
// make_graph / generate / assign_orientation_relaxation first.
//
// Nothing here is new: every population is produced by the SAME call the owning app makes
// (see the table in data/corpus/README.md), and the archived sigma_def orientations are
// read from the SAME directories the apps read them from. The only computation this app
// adds is `checkerboard` (K2c's alternative sigma) and, for the regime population, the
// method::orientation_defect call kill_regime.cpp itself makes because that population has
// no archive.
//
//   freeze_corpus --out <dir> [--only k1a,native200,...] [--recompute-def]
//
// --only          comma-separated subset of the population names below (default: all).
// --recompute-def for yield_fresh: rerun method::orientation_defect and record whether
//                 it reproduces the archived sigma_def (slow, minutes per graph at F~800).
#include <algorithm>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <random>
#include <set>
#include <string>
#include <vector>
#include "kill_common.hpp"
#include "method/design.hpp"

using namespace kiri;
using namespace kiri::kill;
using json = nlohmann::json;

namespace {

json mesh_json_full(const Mesh& m) { return json::parse(mesh_to_json_string(m)); }

// Reads an archived sigma for `m0` exactly as the kill apps do (kind_id.json under dir,
// accepted only when face and vertex counts match). Empty when absent or mismatched.
std::vector<int> archived_sigma(const std::string& dir, const std::string& kind, int id,
                                const Mesh& m0) {
  const std::string p = dir + "/" + kind + "_" + std::to_string(id) + ".json";
  if (!std::filesystem::exists(p)) return {};
  const Mesh sm = load_mesh_json(p);
  if (static_cast<int>(sm.sigma.size()) != m0.n_faces() || sm.X.size() != m0.X.size())
    return {};
  return sm.sigma;
}

// One population row. `mesh.orientation` is sigma_mc (what make_graph leaves in
// mesh.sigma); sigma_mc is repeated as its own field for convenience.
json graph_row(const Graph& g, int min_faces, int max_faces, int n_cap) {
  Mesh m = g.mesh;
  m.build_topology();
  json r;
  r["id"] = g.id;
  r["kind"] = g.kind;
  r["ok"] = g.ok;
  r["make_graph_args"] = {min_faces, max_faces, n_cap};
  r["seed"] = 1000003u * static_cast<unsigned>(g.id) + 20260903u;
  r["N"] = m.n_vertices();
  r["F"] = m.n_faces();
  r["E"] = m.n_edges();
  r["n_interior"] = m.n_interior_vertices();
  r["mesh"] = mesh_json_full(m);
  r["sigma_mc"] = m.sigma;
  r["sigma_checker"] = checkerboard(m);
  r["sigma_def"] = nullptr;
  r["sigma_def_source"] = nullptr;
  return r;
}

void attach_archived_def(json& r, const std::string& dir, const Graph& g) {
  const std::vector<int> sd = archived_sigma(dir, g.kind, g.id, g.mesh);
  if (sd.empty()) return;
  r["sigma_def"] = sd;
  r["sigma_def_source"] = dir + "/" + g.kind + "_" + std::to_string(g.id) + ".json";
}

void dump(const json& j, const std::string& path) {
  write_json(j, path);
  std::cout << "wrote " << path << " (" << std::filesystem::file_size(path) / 1024
            << " KiB)\n";
}

// ids lo..hi (inclusive) through make_graph(id, minf, maxf, ncap); rows for every id,
// including the ones make_graph could not build (ok=false, empty mesh) so the id -> row
// map is total.
json id_range(int lo, int hi, int minf, int maxf, int ncap, const std::string& defdir) {
  json out = json::array();
  for (int id = lo; id <= hi; ++id) {
    const Graph g = make_graph(id, minf, maxf, ncap);
    json r = graph_row(g, minf, maxf, ncap);
    if (g.ok && !defdir.empty()) attach_archived_def(r, defdir, g);
    out.push_back(std::move(r));
    if ((id - lo) % 25 == 24) std::cerr << "  [" << lo << ".." << hi << "] id " << id << "\n";
  }
  return out;
}

// The K5 / K6 / K9 / K9b / K9c / T1 / B4 / native200 / yield(k9c) population: walk ids
// 0..399, keep a graph iff make_graph succeeds AND results/kill/k5/sigma has a matching
// sigma_def, stop after 200. gidx is the population index those apps shard on.
json native200(const std::string& defdir, int n, int maxf) {
  json out = json::array();
  int gidx = -1;
  for (int id = 0; id < 400 && gidx + 1 < n; ++id) {
    const Graph g = make_graph(id, 100, maxf, 1400);
    if (!g.ok) continue;
    const std::vector<int> sd = archived_sigma(defdir, g.kind, id, g.mesh);
    if (sd.empty()) continue;
    ++gidx;
    json r = graph_row(g, 100, maxf, 1400);
    attach_archived_def(r, defdir, g);
    r["gidx"] = gidx;
    out.push_back(std::move(r));
    if (gidx % 25 == 24) std::cerr << "  [native200] gidx " << gidx << " (id " << id << ")\n";
  }
  return out;
}

// kill_regime.cpp: ids 2000..2099, make_graph(id, 20, 100, 1400), sigma_def recomputed
// with method::orientation_defect(m0, sigma_mc, 20 |F|, 7000 + id) (no archive exists).
json regime() {
  json out = json::array();
  for (int id = 2000; id <= 2099; ++id) {
    const Graph g = make_graph(id, 20, 100, 1400);
    json r = graph_row(g, 20, 100, 1400);
    if (g.ok) {
      Mesh m0 = g.mesh;
      m0.build_topology();
      const method::DefectOrientationResult dor =
          method::orientation_defect(m0, m0.sigma, 20 * m0.n_faces(), 7000u + id);
      const bool def_ok = dor.ok && static_cast<int>(dor.sigma.size()) == m0.n_faces();
      // kill_regime falls back to sigma_mc when the search fails; record which happened.
      r["sigma_def"] = def_ok ? dor.sigma : m0.sigma;
      r["sigma_def_source"] = def_ok ? "method::orientation_defect(m0, sigma_mc, 20*F, 7000+id)"
                                     : "fallback: sigma_mc (orientation_defect failed)";
      r["sigma_def_ok"] = def_ok;
      r["sigma_def_defect"] = dor.defect;
    }
    out.push_back(std::move(r));
  }
  return out;
}

// kill_yield.cpp --mode fresh: ids 1000..1099, make_graph(id, 100, 800, 1400), sigma_def
// from results/yield/fresh_sigma (written by that run with orientation_defect, cap 20 |F|,
// seed 7000 + id).
json yield_fresh(const std::string& defdir, bool recompute) {
  json out = json::array();
  for (int id = 1000; id <= 1099; ++id) {
    const Graph g = make_graph(id, 100, 800, 1400);
    json r = graph_row(g, 100, 800, 1400);
    if (g.ok) {
      attach_archived_def(r, defdir, g);
      if (recompute) {
        Mesh m0 = g.mesh;
        m0.build_topology();
        Timer t;
        const method::DefectOrientationResult dor =
            method::orientation_defect(m0, m0.sigma, 20 * m0.n_faces(), 7000u + id);
        r["sigma_def_recomputed"] = dor.sigma;
        r["sigma_def_recompute_ok"] = dor.ok;
        r["sigma_def_recompute_matches"] =
            !r["sigma_def"].is_null() && r["sigma_def"].get<std::vector<int>>() == dor.sigma;
        std::cerr << "  [yield_fresh] id " << id << " recompute " << (dor.ok ? "ok" : "FAIL")
                  << " match=" << r["sigma_def_recompute_matches"] << " " << t.s() << " s\n";
      }
    }
    out.push_back(std::move(r));
    if (id % 25 == 24) std::cerr << "  [yield_fresh] id " << id << "\n";
  }
  return out;
}

// kill_k2b.cpp: ids 0..399 through make_graph(id, 100, 160, 1400), first 8 that build.
json k2b() {
  json out = json::array();
  int taken = 0;
  for (int id = 0; id < 400 && taken < 8; ++id) {
    const Graph g = make_graph(id, 100, 160, 1400);
    if (!g.ok) continue;
    json r = graph_row(g, 100, 160, 1400);
    r["gidx"] = taken++;
    out.push_back(std::move(r));
  }
  return out;
}

// tests/derivation_tests.cpp l1_corpus(): make_graph(id, 18, 46, 220) for id = 0..79,
// keeping every graph that builds (the test stops once it has 60 cases in total).
// derivations/scratch/check_l1 (`check_l1 140 20`) walks the same call for ids 0..139.
json derivation_l1(int n_ids = 80) {
  json out = json::array();
  for (int id = 0; id < n_ids; ++id) {
    const Graph g = make_graph(id, 18, 46, 220);
    if (!g.ok) continue;
    out.push_back(graph_row(g, 18, 46, 220));
  }
  return out;
}

// kill_scaling.cpp build_graph(kind, sites, seed, true) on the grid of
// code/scripts/run_scaling.sh: targets {50,...,5000}, delaunay sites = target/2, voronoi
// sites = target, seeds 0..2.
json scaling() {
  json out = json::array();
  const std::vector<int> targets{50, 100, 200, 500, 1000, 2000, 5000};
  for (const std::string kind : {"delaunay", "voronoi"}) {
    for (int t : targets) {
      const int sites = kind == "delaunay" ? t / 2 : t;
      for (int seed = 0; seed < 3; ++seed) {
        const unsigned s = 1000003u * static_cast<unsigned>(seed) + 20260908u +
                           7919u * static_cast<unsigned>(sites);
        std::mt19937 rng(s);
        Mesh m = generate(kind, {static_cast<double>(sites), 40.0}, rng);
        m.build_topology();
        m.sigma = assign_orientation_relaxation(m, rng, 4, 300, 90).sigma;
        json r;
        r["kind"] = kind;
        r["target_faces"] = t;
        r["sites"] = sites;
        r["seed"] = seed;
        r["rng_seed"] = s;
        r["N"] = m.n_vertices();
        r["F"] = m.n_faces();
        r["E"] = m.n_edges();
        r["mesh"] = mesh_json_full(m);
        r["sigma_mc"] = m.sigma;
        out.push_back(std::move(r));
        std::cerr << "  [scaling] " << kind << " sites=" << sites << " seed=" << seed
                  << " F=" << m.n_faces() << "\n";
      }
    }
  }
  return out;
}

// kill_common::reference_cases(): the eight Phase-2 cases with their sigma.
json reference_cases_json() {
  json out = json::array();
  for (const RefCase& c : reference_cases()) {
    Mesh m = c.mesh;
    m.build_topology();
    json r;
    r["name"] = c.name;
    r["periodic"] = c.periodic;
    r["N"] = m.n_vertices();
    r["F"] = m.n_faces();
    r["mesh"] = mesh_json_full(m);
    r["sigma"] = m.sigma;
    out.push_back(std::move(r));
  }
  return out;
}

// kill_common::deployable_population(samples, radius): (8, 0.2) is the default used by
// kill_k1a / kill_k2b, (2, 0.2) by kill_k8a, (20, 0.2) and (20, 0.35) by kill_e1 part (A).
json deployable_json(int samples = 8, double radius = 0.2) {
  json out = json::array();
  for (const DeployableConfig& d : deployable_population(samples, radius)) {
    json r;
    r["name"] = d.name;
    r["family"] = d.family;
    r["dim_null"] = d.dim_null;
    r["sample"] = d.sample;
    r["mesh"] = mesh_json_full(d.mesh);
    r["X"] = json::array();
    for (const auto& x : d.X) r["X"].push_back({x.x(), x.y()});
    out.push_back(std::move(r));
  }
  return out;
}

std::string g17(double v) {
  char b[64];
  std::snprintf(b, sizeof b, "%.17g", v);
  return b;
}

// Test vectors for a bit-exact std::mt19937 + libc++ distribution port.
json mt_vectors() {
  json out;
  out["__cplusplus"] = static_cast<long long>(__cplusplus);
#ifdef _LIBCPP_VERSION
  out["_LIBCPP_VERSION"] = static_cast<long long>(_LIBCPP_VERSION);
#else
  out["_LIBCPP_VERSION"] = nullptr;
#endif
#ifdef __clang_version__
  out["compiler"] = std::string("clang ") + __clang_version__;
#elif defined(__VERSION__)
  out["compiler"] = __VERSION__;
#endif
  out["note"] =
      "Each list is drawn from a FRESH std::mt19937(seed); the distribution objects are "
      "fresh too. real draws are printed with %.17g and also stored as JSON numbers. "
      "shuffle_0_19 is std::shuffle on [0..19] with a fresh engine.";
  out["seeds"] = json::array();
  const std::vector<unsigned> seeds{0u, 1u, 5489u, 20260903u, 1000003u * 7u + 20260903u};
  for (unsigned seed : seeds) {
    json s;
    s["seed"] = seed;
    {
      std::mt19937 e(seed);
      json v = json::array();
      for (int i = 0; i < 20; ++i) v.push_back(static_cast<unsigned long long>(e()));
      s["raw_u32"] = v;
    }
    auto reals = [&](double lo, double hi, const char* key) {
      std::mt19937 e(seed);
      std::uniform_real_distribution<double> U(lo, hi);
      json v = json::array(), vs = json::array();
      for (int i = 0; i < 20; ++i) {
        const double d = U(e);
        v.push_back(d);
        vs.push_back(g17(d));
      }
      s[key] = v;
      s[std::string(key) + "_str"] = vs;
    };
    reals(0.0, 1.0, "uniform_real_0_1");
    reals(0.0, 40.0, "uniform_real_0_40");
    reals(0.0, 2 * M_PI, "uniform_real_0_2pi");
    {
      std::mt19937 e(seed);
      std::uniform_int_distribution<int> D(0, 99);
      json v = json::array();
      for (int i = 0; i < 20; ++i) v.push_back(D(e));
      s["uniform_int_0_99"] = v;
    }
    {
      std::mt19937 e(seed);
      std::uniform_int_distribution<int> D(0, 1);
      json v = json::array();
      for (int i = 0; i < 20; ++i) v.push_back(D(e));
      s["uniform_int_0_1"] = v;
    }
    {
      std::mt19937 e(seed);
      std::normal_distribution<double> G(0.0, 1.0);
      json v = json::array(), vs = json::array();
      for (int i = 0; i < 20; ++i) {
        const double d = G(e);
        v.push_back(d);
        vs.push_back(g17(d));
      }
      s["normal_0_1"] = v;
      s["normal_0_1_str"] = vs;
    }
    {
      std::mt19937 e(seed);
      std::vector<int> p(20);
      for (int i = 0; i < 20; ++i) p[i] = i;
      std::shuffle(p.begin(), p.end(), e);
      s["shuffle_0_19"] = p;
    }
    out["seeds"].push_back(std::move(s));
  }
  return out;
}

}  // namespace

int main(int argc, char** argv) {
  std::string out = "/Users/emredayangac/Documents/kirigami-julia/data/corpus";
  std::string k5_sigma = "results/kill/k5/sigma", fresh_sigma = "results/yield/fresh_sigma";
  std::set<std::string> only;
  bool recompute_def = false;
  for (int i = 1; i < argc; ++i) {
    const std::string a = argv[i];
    if (a == "--out" && i + 1 < argc) out = argv[++i];
    else if (a == "--k5-sigma" && i + 1 < argc) k5_sigma = argv[++i];
    else if (a == "--fresh-sigma" && i + 1 < argc) fresh_sigma = argv[++i];
    else if (a == "--recompute-def") recompute_def = true;
    else if (a == "--only" && i + 1 < argc) {
      std::stringstream ss(argv[++i]);
      std::string tok;
      while (std::getline(ss, tok, ',')) only.insert(tok);
    }
  }
  auto want = [&](const std::string& k) { return only.empty() || only.count(k); };
  std::filesystem::create_directories(out);
  Timer wall;

  if (want("mt19937")) dump(mt_vectors(), out + "/mt19937_vectors.json");
  if (want("reference_cases")) dump(reference_cases_json(), out + "/reference_cases_8.json");
  if (want("derivation_l1")) dump(derivation_l1(80), out + "/derivation_l1_small.json");
  if (want("derivation_l1_140")) dump(derivation_l1(140), out + "/derivation_l1_140.json");
  if (want("k2b")) dump(k2b(), out + "/k2b_random_8.json");
  if (want("regime")) dump(regime(), out + "/regime_100.json");
  if (want("k1a")) dump(id_range(0, 199, 100, 800, 1400, k5_sigma), out + "/k1a_200.json");
  if (want("native200")) dump(native200(k5_sigma, 200, 800), out + "/native200.json");
  if (want("yield_fresh"))
    dump(yield_fresh(fresh_sigma, recompute_def), out + "/yield_fresh_100.json");
  if (want("e1")) dump(id_range(0, 899, 100, 800, 1400, k5_sigma), out + "/e1_900.json");
  if (want("deployable")) dump(deployable_json(8, 0.2), out + "/deployable_population.json");
  if (want("deployable_variants")) {
    dump(deployable_json(2, 0.2), out + "/deployable_population_2_0.2.json");
    dump(deployable_json(20, 0.2), out + "/deployable_population_20_0.2.json");
    dump(deployable_json(20, 0.35), out + "/deployable_population_20_0.35.json");
  }
  if (want("scaling")) dump(scaling(), out + "/scaling_42.json");
  if (want("k2c")) dump(id_range(0, 499, 100, 5000, 1400, ""), out + "/k2c_500.json");
  if (want("k3a")) dump(id_range(0, 499, 100, 5000, 1 << 30, ""), out + "/k3a_500.json");

  std::cout << "done in " << wall.s() << " s\n";
  return 0;
}
