#include <iomanip>
#include <iostream>
#include <fstream>

using namespace std;

int
main(int argc, char* argv[])
{
  if (argc < 3) {
    std::cout << argv[0] << " error: provide at least two command line arguments\n";
    return 1;
  }


  ofstream ofs(argv[argc - 1], ofstream::out);

  for (int n = 1; n < argc - 1; n++) {
    ifstream ifs(argv[n], ifstream::in);

    char c = ifs.get();

    while (ifs.good()) {
      ofs << c;
      c = ifs.get();
    }
    ifs.close();
  }

  ofs.close();

  return 0;
}
