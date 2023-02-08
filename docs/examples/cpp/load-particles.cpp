#include <hippio.h>

using namespace HIPP;
using namespace HIPP::IO;

int main(int argc, char const *argv[]) {
    H5File f("simulation.hdf5", "r");

    H5Group header = f.open_group("Header");

    double h, Omg0;
    header.open_attr("HubbleParam").read(&h);
    header.open_attr("Omega0").read(&Omg0);
    pout << "hubble parameter = ", h, ", matter density = ", Omg0, endl;


    vector<double> zs;
    header.open_dataset("Redshifts").read(zs);
    pout << "Redshifts = {", zs, "}", endl;


    auto parts = f.open_group("Snapshots/76/0/PartType1");
    vector<string> keys = parts.keys();
    pout << "Names of fields: {", keys, "}\n";
    
    long long n_parts;
    parts.open_attr("NumPart_ThisFile").read(&n_parts);
    pout << "No. of particles = ", n_parts, endl;

    vector<std::array<float, 3>> xs(n_parts);
    parts.open_dataset("Coordinates").read(&xs[0][0]);
    pout << "Positions of the first five particles: \n";
    for(int i=0; i<5; ++i)
        pout << "  {", xs[i], "}\n";

    return 0;
}

