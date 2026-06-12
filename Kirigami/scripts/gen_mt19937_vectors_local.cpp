#include <random>
#include <algorithm>
#include <cstdio>
#include <cstdint>
#include <cstring>
#include <vector>
// Emits reference vectors for the Julia port of std::mt19937 + libc++ distributions.
static void hex(double x){ uint64_t b; std::memcpy(&b,&x,8); std::printf("\"0x%016llx\"",(unsigned long long)b); }
int main(){
  const unsigned seeds[] = {0u,1u,5489u,20260903u};
  const int N = 20;
  std::printf("{\n\"libcpp_version\": %d,\n", _LIBCPP_VERSION);
  std::printf("\"textbook_10000th\": %u,\n", [](){ std::mt19937 g; for(int i=0;i<9999;++i) g(); return g(); }());
  std::printf("\"seeds\": [\n");
  for(size_t s=0;s<4;++s){
    unsigned seed=seeds[s];
    std::printf("{\"seed\": %u,\n", seed);
    { std::mt19937 g(seed); std::printf("\"u32\": ["); for(int i=0;i<N;++i) std::printf("%s%u", i?",":"", g()); std::printf("],\n"); }
    { std::mt19937 g(seed); std::printf("\"canonical\": ["); for(int i=0;i<N;++i){ if(i) std::printf(","); hex(std::generate_canonical<double,53>(g)); } std::printf("],\n"); }
    { std::mt19937 g(seed); std::uniform_real_distribution<double> U(-1.5,1.5); std::printf("\"uniform_real\": {\"a\": -1.5, \"b\": 1.5, \"values\": ["); for(int i=0;i<N;++i){ if(i) std::printf(","); hex(U(g)); } std::printf("]},\n"); }
    { std::mt19937 g(seed); std::uniform_real_distribution<double> U(0.05,1.4); std::printf("\"uniform_real2\": {\"a\": 0.05, \"b\": 1.4, \"values\": ["); for(int i=0;i<N;++i){ if(i) std::printf(","); hex(U(g)); } std::printf("]},\n"); }
    // several ranges: power-of-two-minus-one, small, wide, single
    const int ranges[][2] = {{0,6},{0,7},{0,99},{-3,3},{5,5},{0,1000000}};
    std::printf("\"uniform_int\": [");
    for(int r=0;r<6;++r){ std::mt19937 g(seed); std::uniform_int_distribution<int> D(ranges[r][0],ranges[r][1]);
      std::printf("%s{\"a\": %d, \"b\": %d, \"values\": [", r?",":"", ranges[r][0], ranges[r][1]);
      for(int i=0;i<N;++i) std::printf("%s%d", i?",":"", D(g)); std::printf("]}"); }
    std::printf("],\n");
    { std::mt19937 g(seed); std::vector<int> v(10); for(int i=0;i<10;++i) v[i]=i; std::shuffle(v.begin(),v.end(),g);
      std::printf("\"shuffle\": {\"n\": 10, \"perm\": ["); for(int i=0;i<10;++i) std::printf("%s%d", i?",":"", v[i]); std::printf("], \"next_u32\": %u}\n", g()); }
    std::printf("}%s\n", s<3?",":"");
  }
  std::printf("]\n}\n");
}
