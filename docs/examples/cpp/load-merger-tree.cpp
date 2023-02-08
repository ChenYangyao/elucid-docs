#include <hippio.h>

using namespace HIPP;
using namespace HIPP::IO;


struct Subhalo {
    float mass;
    float x[3];
    int fpro_id, cent_id;
    int snap;
};

int main(int argc, char const *argv[]) {
    H5File f("simulation.hdf5", "r");
    auto g = f.open_group("Trees/SubLink/0");
    
    vector<Subhalo> forest = H5XTable<Subhalo>(
        "M_TopHat", &Subhalo::mass,                 "Pos", &Subhalo::x,
        "FirstProgenitor", &Subhalo::fpro_id,       "FirstHaloInFOFgroup", &Subhalo::cent_id,
        "SnapNum", &Subhalo::snap
    ).read(g.open_group("Subhalos"));

    vector<int> offs, lens;
    g.open_dataset("Header/HaloOffsetInTree").read(offs);
    g.open_dataset("Header/NumHalosInTree").read(lens);

    int nhalos = forest.size(), ntrees = offs.size();
    pout << "No. of subhalos/trees in this chunk = ", nhalos, '/', ntrees, 
        endl;

    int off = offs[0], len = lens[0];
    const auto *tree = &forest[off];

    int root_id = -1;
    for(int i=0; i<len; ++i){
        auto &h = tree[i];
        if( h.snap == 100 && h.cent_id == i && h.mass > 1.0e2 ) {
            root_id = i; break;
        }
    }
    pout << "A root index at z=0 : ", root_id, endl;

    int pro_id = root_id;
    while( pro_id >= 0 ) {
        auto &pro = tree[pro_id];
        if( pro.mass < tree[root_id].mass * 0.5 ) {
            pout << "Half mass formation snapshot: ", pro.snap, endl;
            break;
        }
        pro_id = pro.fpro_id;
    }

    return 0;
}

