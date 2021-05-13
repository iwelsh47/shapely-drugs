#include <Eigen/Dense>
#include <iostream>
#include <iomanip>
#include <Shapely/coordinates.hpp>
#include <Shapely/molecule.hpp>
#include <vector>
#include <bitset>
#include <array>
#include <sul/dynamic_bitset.hpp>
#include <Shapely/serialisation.hpp>
#include <Shapely/compressed_bitset.hpp>
#include <Shapely/utils.hpp>
#include <filesystem>
#include <fstream>

#if 0
#include <GraphMol/GraphMol.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/FileParsers/MolWriters.h>
#include <GraphMol/Descriptors/MolDescriptors.h>
#include <GraphMol/MolAlign/AlignMolecules.h>
#include <GraphMol/DistGeomHelpers/Embedder.h>
#include <GraphMol/ForceFieldHelpers/MMFF/MMFF.h>
#include <GraphMol/ForceFieldHelpers/UFF/UFF.h>
#endif

#include <Shapely/FFTProcess.hpp>
#include <Shapely/grid.hpp>

namespace fs = std::filesystem;

CEREAL_CLASS_VERSION(shapely::Grid, 1);

void print_header(std::ofstream& os) {
  os << std::setprecision(6) << std::scientific;
  os << std::setw(12) << "Molecule A" << " "
     << std::setw(12) << "Molecule B" << " "
     << std::setw(3)  << "STY" << " "
     << std::setw(10) << "Volume A" << " "
     << std::setw(10) << "Volume B" << " "
     << std::setw(10) << "A & B" << " "
     << std::setw(10) << "A | B" << " "
     << std::setw(10) << "A - B" << " "
     << std::setw(10) << "B - A" << " "
     << std::setw(15) << "ESP (outer)" << " "
     << std::setw(15) << "ESP (A only)" << " "
     << std::setw(15) << "ESP (B only)" << '\n';
}

void print_line(std::ofstream& os, stdx::string a, stdx::string b, const shapely::Overlap& best_overlap, uint8_t best_style) {
  os << std::setw(12) << a << " "
     << std::setw(12) << b << " "
     << std::setw(3)  << (int32_t)best_style << " "
     << std::setw(10) << best_overlap.VA << " "
     << std::setw(10) << best_overlap.VB << " "
     << std::setw(10) << best_overlap.AandB << " "
     << std::setw(10) << best_overlap.AplusB << " "
     << std::setw(10) << best_overlap.AminusB << " "
     << std::setw(10) << best_overlap.BminusA << " "
     << std::setw(15) << best_overlap.esp_outerV << " "
     << std::setw(15) << best_overlap.esp_Va << " "
     << std::setw(15) << best_overlap.esp_Vb << std::endl;
}

void run_compare(std::ofstream& os, const shapely::Grid& g1, const shapely::Grid& g2) {
  
  // Get the best overlap
//  auto best = g1.CalculateBestOverlap(g2);
  // Get all the overlaps
  for (uint8_t style = 0; style < g1.volume.size(); ++style) {
    shapely::Overlap overlap = g1.CalculateOverlap(g2, style);
    print_line(os, g1.name, g2.name, overlap, style);
  }
  os.flush();
  
  // Print out the best data
//  print_line(os, g1.name, g2.name, best.first, best.second);
}

void LoadFile(const fs::path& file, std::vector<shapely::Grid>& data) {
  shapely::Molecule discard;
  using Deserialiser = cereal::PortableBinaryInputArchive;
  std::ifstream in(file, std::ios::in | std::ios::binary);
  Deserialiser loader(in);
  decltype(data.size()) num_molecules;
  loader(num_molecules);
  data.resize(num_molecules);
  for (uint64_t i = 0; i < num_molecules; ++i) {
    loader(discard, data[i]);
  }
}

int file_check_main(int argc, const char* argv[]) {
  fs::path dump_dir_root("/home/ivan/shapely-drugs/data/SHAPE_DUMPS/");
  for (fs::directory_entry file : fs::directory_iterator(dump_dir_root)) {
    if (file.path().extension() != ".dump") { continue; }
    std::cout << "Loading " << file.path() << "...";
    std::vector<shapely::Grid> data;
    try {
	LoadFile(file.path(), data);
    } catch (cereal::Exception) {
	    std::cout << " FAIL." << std::endl;
	    continue;
    }
    std::cout << " SUCCESS!." << std::endl;
  }
  return 0;
}

int scoring_main(int argc, const char* argv[]) {
#if 1
  fs::path dump_dir_root("/Users/ivanwelsh/GitHub/shapely-drugs/data/ESR_DRUGS/");
  fs::path output_dir_root("/Users/ivanwelsh/GitHub/shapely-drugs/data/analysis/ESR_DRUGS/");
#else
  fs::path dump_dir_root("/home/ivan/shapely-drugs/data/ESR_DRUGS/");
  fs::path output_dir_root("/home/ivan/shapely-drugs/data/analysis/ESR_DRUGS/");
#endif
//  std::vector<stdx::string> targets{
//    // Desired targets
//    "er_agonist", "er_antagonist", //"ache", "comt", "ar", "gr", "mr", "pde5", "pnp", "pr",
//    // Extra targets
//    //"ada", "alr2", "ampc",  "cdk2",  "cox1", "cox2", "dhfr", "egfr", "ace", "fgfr1", "fxa", "gart",
//    //"gpb", "hivpr", "hivrt", "hmga", "hsp90", "inha", "na", "p38", "parp", "pdgfrb", "ppar_gamma",
//    //"rxr_alpha", "sahh", "src", "thrombin", "tk", "trypsin", "vegfr2"
//  };

//  for (stdx::string& target : targets) {
//    if (target != "comt") { continue; }
    fs::path ligands_file = dump_dir_root / "all_drugs.dump";
//    fs::path decoys_file = dump_dir_root / (target + "_decoys.dump");
    
    // Make sure we have somewhere to save the output
    fs::path out_root = output_dir_root /  "scores/";
    fs::create_directories(out_root);
    
    // Load the target ligand + decoy data
    std::vector<shapely::Grid> ligands, decoys;
    LoadFile(ligands_file, ligands);
 //   LoadFile(decoys_file, decoys);
    
    // Parallelise on first iteration through the ligands
#pragma omp parallel for schedule(static, 1)
    for (uint64_t ligand_idx = 0; ligand_idx < ligands.size(); ++ligand_idx) {
      shapely::Grid* grid_i = new shapely::Grid(ligands[ligand_idx]);
      fs::path outfile = out_root / (grid_i->name + ".txt");
      if (fs::exists(outfile) && fs::file_size(outfile) > 1048576) {
        delete grid_i;
        continue;
      }
      
      grid_i->Expand();
      // Open the ligand-ligand output file
#if defined(_OPENMP)
#pragma omp critical
      std::cout << "Running ligand compare with " << grid_i->name << " on thread " << omp_get_thread_num() << std::endl;
#else
      std::cout << "Running ligand compare with " << grid_i->name << std::endl;
#endif
      std::ofstream ligand_out(outfile);
      print_header(ligand_out);
      
      // Iterate through the ligands, doing ligand-ligand comparisons
      for (const shapely::Grid& grid_raw : ligands) {
        if (grid_raw.name == grid_i->name) { continue; }
        shapely::Grid* grid_j = new shapely::Grid(grid_raw);
        grid_j->Expand();
        run_compare(ligand_out, *grid_i, *grid_j);
        delete grid_j;
      }
      ligand_out.close();
      
      // Iterate through the decoys, doing decoy-ligand comparisons
      // Open the ligand-decoy output file
//#if defined(_OPENMP)
//#pragma omp critical
//      std::cout << "Running decoy compare with " << target << "/" << grid_i->name << " on thread " << omp_get_thread_num() << std::endl;
//#else
//      std::cout << "Running decoy compare with " << target << "/" << grid_i->name << std::endl;
//#endif
//      std::ofstream decoy_out(out_root / (grid_i->name + "_decoys.txt"));
//      print_header(decoy_out);
//      
//      // Iterate through the decoys, doing ligand-decoy comparisons
//      for (const shapely::Grid& grid_raw : decoys) {
//        shapely::Grid* grid_j = new shapely::Grid(grid_raw);
//        grid_j->Expand();
//        run_compare(decoy_out, *grid_i, *grid_j);
//        delete grid_j;
//      }
//      decoy_out.close();
//      delete grid_i;
    }

//  }
  return 0;
}

int prep_main(int argc, const char * argv[]) {
  using Serialiser = cereal::PortableBinaryOutputArchive;
  using Deserialiser = cereal::PortableBinaryInputArchive;
#if 0
  fs::path atoms_dir("/Users/ivanwelsh/GitHub/shapely-drugs/data/ISA_ATOMS/");
  fs::path dump_dir_root("/Volumes/RealData/SHAPE_DUMPS/");
#else
  fs::path dump_dir_root("/home/ivan/shapely-drugs/data/ESR_DRUGS/");
  fs::path atoms_dir("/home/ivan/shapely-drugs/data/ESR_DRUGS/ISA_ATOMS/");
#endif

 // for (const fs::directory_entry& file : fs::directory_iterator(atoms_dir)) {
 //   if (!file.is_directory()) { continue;; }
 //   fs::path isaq_dir = file.path();
 //   fs::path dump_path = dump_dir_root / isaq_dir.filename();
 //   dump_path.replace_extension("dump");

 //   std::cout << "Prepping " << dump_path.stem() << std::endl;
    std::vector<fs::path> files;
    fs::path dump_path = dump_dir_root / "all_drugs.dump";
    for (const fs::directory_entry& isaq_file : fs::directory_iterator(atoms_dir)) {
      if (!isaq_file.is_regular_file()) { continue; }
      if (isaq_file.path().extension() != ".isaq") { continue; }
      files.emplace_back(isaq_file.path());
    }
    
    std::ofstream dumpstream(dump_path, std::ios::out | std::ios::binary);
    Serialiser dumper(dumpstream);
    dumper(files.size());
    
#pragma omp parallel for schedule(guided, 1)
    for (uint64_t i = 0; i < files.size(); ++i) {
      fs::path infile = files[i];
      shapely::Molecule* mol = new shapely::Molecule;
      mol->LoadFile(infile);
      shapely::Grid* grid = new shapely::Grid(*mol);
      grid->Compress();
#pragma omp critical
      dumper(*mol, *grid);
      delete mol;
      delete grid;
    }
 // }
  
  return 0;
}

// Loads all the dump files in SHAPE_DUMPS and resaves with new save version.
int update_savedfiles(int argc, const char* argv[]) {
  using Deserialiser = cereal::BinaryInputArchive;
  using Serialiser = cereal::BinaryOutputArchive;
  
  fs::path dump_dir_root("/Volumes/RealData/SHAPE_DUMPS/");
  for (const fs::directory_entry& file : fs::recursive_directory_iterator(dump_dir_root)) {
    if (file.path().extension() != ".dump") { continue; }
    fs::path new_name = file.path();
    new_name.replace_extension(".dump2");
    fs::rename(file.path(), new_name);
    
    std::cout << std::setw(25) << new_name.stem();
    
    std::vector<fs::path>::size_type num_systems;
    std::ifstream in(new_name, std::ios::in | std::ios::binary);
    std::ofstream out(file.path(), std::ios::out | std::ios::binary);
    {
      Deserialiser loader(in);
      loader(num_systems);
      
      Serialiser dumper(out);
      dumper(num_systems);
      // Load the file
      shapely::Molecule mol_i; shapely::Grid grid_i;
      for (uint64_t i = 0; i < num_systems; ++i) {
        loader(mol_i, grid_i);
        dumper(mol_i, grid_i);
      }
    }
    std::cout << ": Initial = " << fs::file_size(new_name) << ", Final = " << fs::file_size(file.path()) << '\n';
    fs::remove(new_name);
  }
  
  return 0;
}

int print_coordinates(int argc, const char* argv[]) {
  using Deserialiser = cereal::BinaryInputArchive;
  fs::path dump_dir_root("/Volumes/RealData/SHAPE_DUMPS/");
  fs::path output_dir_root("/Users/ivanwelsh/GitHub/shapely-drugs/data/analysis/");
  std::vector<stdx::string> targets{
    "er_agonist", "er_antagonist",
    "ace", "ache", "ada", "alr2", "ampc", "ar", "cdk2", "comt", "cox1", "cox2", "dhfr", "egfr",
    "fgfr1", "fxa", "gart", "gpb", "gr", "hivpr", "hivrt", "hmga",
    "hsp90", "inha", "mr", "na", "p38", "parp", "pde5", "pdgfrb", "pnp", "ppar_gamma", "pr",
    "rxr_alpha", "sahh", "src", "thrombin", "tk", "trypsin", "vegfr2"
  };
  
  for (const stdx::string& target : targets) {
    std::cout << "Working on " << target << std::endl;
    fs::path ligands = dump_dir_root / (target + "_ligands.dump");
    shapely::Molecule mol; shapely::Grid grid;
    
    std::ifstream in(ligands, std::ios::in | std::ios::binary);
    Deserialiser data(in);
    // Double iterate through the ligand files
    while (true) {
      // Load the ligand data, one by one
      try {
        data(mol, grid);
      } catch (cereal::Exception) {
        break;
      }
//      if (mol.name != "ZINC00001219") { continue; }
      
      // Make sure we have somewhere to save the output
      fs::path out_root = output_dir_root / (target + "/ligand_structures");
      fs::create_directories(out_root);
      
      std::ofstream out(out_root/(mol.name + ".xyz"));
      mol.PrintXYZ(out);
      out.close();
      out.open(out_root/("vol_" + mol.name + ".xyz"));
      grid.PrintXYZ(out, shapely::Limits(mol.atoms, mol.radii.data()), 0);
    }
  }
  return 0;
}

int test_main(int argc, const char* argv[]) {
  using namespace shapely;
  fs::path testfile("/Users/ivanwelsh/GitHub/shapely-drugs/data/test.isaq");
  Molecule mol;
  mol.LoadFile(testfile);
  
  std::ofstream out_mol("/Users/ivanwelsh/GitHub/shapely-drugs/data/test.xyz");
  mol.PrintXYZ(out_mol);
  std::ofstream out_vol("/Users/ivanwelsh/GitHub/shapely-drugs/data/test_vol2.xyz");
  
  Grid grid(mol);
  grid.PrintXYZ(out_vol, Limits(mol.atoms, mol.radii.data()), 0);
  
  return 0;
}

#if 0
int main_rdkittest(int argc, const char* argv[]) {
  
  fs::path sdf_root("/Users/ivanwelsh/GitHub/shapely-drugs/data/DUD/");
  fs::path isaq_root("/Users/ivanwelsh/GitHub/shapely-drugs/data/ISA_ATOMS/");
  
  std::unique_ptr<RDKit::ROMol> mol;
  
  for (const fs::directory_entry& item : fs::directory_iterator(sdf_root)) {
    // Check for sdf file
    if (item.path().extension() != ".sdf") { continue; }
    fs::path sdf_file = item.path();
    if (sdf_file.stem() != "comt_ligands") { continue; }
    
    // Load the SDF file and iterate through the molecules
    RDKit::SDMolSupplier mol_supplier(sdf_file.native(), true, false, true);
    while (!mol_supplier.atEnd()) {
      mol.reset(mol_supplier.next());
      stdx::string name = mol->getProp<std::string>("_Name");
      fs::path isaq_file = isaq_root / sdf_file.stem() / (name + ".isaq");
      if (!fs::exists(isaq_file)) { continue; }
      uint32_t num_rot_bonds = RDKit::Descriptors::calcNumRotatableBonds(*mol, false);
      std::cout << "Loaded " << name << " which has " << num_rot_bonds << " rotatable bonds." << std::endl;
      
      RDKit::DGeomHelpers::EmbedParameters params(RDKit::DGeomHelpers::ETKDGv3);
      params.useSmallRingTorsions = true;
      params.randomSeed = 123456;
      
      RDKit::INT_VECT mol_conformer_ids = RDKit::DGeomHelpers::EmbedMultipleConfs(*mol, num_rot_bonds * 100, params);
      std::cout << "Generated " << mol_conformer_ids.size() << " conformers." << std::endl;
      
      std::vector<double> rms_list;
      RDKit::MolAlign::alignMolConformers(*mol, nullptr, nullptr, nullptr, false, 50, &rms_list);
      std::sort(rms_list.begin(), rms_list.end());
      for (double& v : rms_list) {
        continue;
      }
    }
    
  }
  
  return 0;
}
#endif

int main_esgrid_test(int argc, const char* argv[]) {
  fs::path root_path("/Users/ivanwelsh/GitHub/shapely-drugs/data/ESR_DRUGS/ISA_ATOMS");
  fs::path target_path = root_path / "Methamphetamine.density";
  std::vector<fs::path> test_paths {
    root_path / "2-Chloromphethamine.density",
    root_path / "2-Fluroamphethamine.density",
    root_path / "2-Fluromethamphethamine.density",
    root_path / "2C-B_from-amphethamine.density"
  };
  
  shapely::Molecule target;
  target.LoadBinaryFile(target_path);
  
  return 0;
}

int main_fftw_test(int argc, const char* argv[]) {
  fs::path root_path("/Users/ivanwelsh/GitHub/shapely-drugs/data/ESR_DRUGS/ISA_ATOMS");
  fs::path target_path = root_path / "Methamphetamine.density";
  std::vector<fs::path> test_paths {
    root_path / "2-Chloromphethamine.density",
    root_path / "2-Fluroamphethamine.density",
    root_path / "2-Fluromethamphethamine.density",
    root_path / "2C-B_from-amphethamine.density",
    root_path / "Methamphetamine.density"
  };
  
  
//  for (fs::directory_entry file : fs::directory_iterator(root_path)) {
//    if (!fs::is_regular_file(file) || file.path().extension() != ".isaq") { continue; }
//    test_paths.push_back(file.path());
//  }
  
  std::map<uint32_t, std::vector<shapely::Score>> results;
  
  shapely::FFTProcess runner(target_path, test_paths);
  runner.Initialise();
  runner.Run(25, results);
  fs::path target_out = target_path; target_out.replace_extension(".target.xyz");
  {
    std::ofstream os(target_out);
    runner.PrintTargetXYZ(os);
  }
  
  for (auto& iv : results) {
    std::cout << "Results of " << test_paths[iv.first].filename() << ":\n";
    for (shapely::Score& score : iv.second) {
      std::cout << std::setw(4) << score.GetRotationType() << ": "
                << std::setw(10) << score.GetRawIndex() << " ("
                << std::setw(3) << score.GetXIndex() << ", "
                << std::setw(3) << score.GetYIndex() << ", "
                << std::setw(3) << score.GetZIndex() << ")  "
                << std::setprecision(6) << std::fixed << score.GetScore() << std::endl;
      {
        std::stringstream ss;
        ss << ".test." << std::setw(4) << std::setfill('0') << score.GetRotationType() << ".xyz";
        fs::path test_out = test_paths[iv.first];
        test_out.replace_extension(ss.str());
        std::ofstream os(test_out);
        runner.PrintScoredTestXYZ(os, iv.first, score);
        if (score.GetRotationType() == 0) {
          ss.str("");
          ss << ".test_initial.xyz";
          os.close();
          test_out = test_paths[iv.first];
          test_out.replace_extension(ss.str());
          os.open(test_out);
          runner.PrintTestXYZ(os, iv.first);
        }
      }
    }
  }
  
  return 0;
}

int main_fftw_allbyall(int argc, const char* argv[]) {
  std::ofstream mat_out("/Users/ivanwelsh/GitHub/shapely-drugs/data/ESR_DRUGS/10_by_10.mat");
  std::ofstream raw_out("/Users/ivanwelsh/GitHub/shapely-drugs/data/ESR_DRUGS/10_by_10.dat");
  mat_out.precision(6);
  
  fs::path root_path("/Users/ivanwelsh/GitHub/shapely-drugs/data/ESR_DRUGS/ISA_ATOMS2");
  std::vector<fs::path> test_files;
  for (const fs::directory_entry& possible : fs::directory_iterator(root_path)) {
    if (possible.path().extension() != ".density") { continue; }
    stdx::string name = possible.path().stem().native();
    if (name == "10_41_Lysergic_acid") { continue; }
    stdx::string start(name.begin(), name.begin() + 3);
    if (start != "10_") { continue; }
    test_files.emplace_back(possible.path());
  }
  std::sort(test_files.begin(), test_files.end());
  
//#pragma omp parallel for
  for (uint64_t i = 0; i < test_files.size(); ++i) {
    fs::path p1 = test_files[i];
    std::cout << "Running target = " << p1.stem() << '\n';
    raw_out << "Results of " << p1.stem() << ":" << std::endl;
    std::map<uint32_t, std::vector<shapely::Score>> results;
    shapely::FFTProcess runner(p1, test_files);
    runner.Initialise();
    runner.Run(1, results);
    
    std::stringstream ss;
    ss << "/Users/ivanwelsh/GitHub/shapely-drugs/data/ESR_DRUGS/TestOverlays/" << p1.stem().native() << ".target.xyz";
    std::ofstream target_out(ss.str());
    runner.PrintTargetXYZ(target_out);
    target_out.close();
    
    for (auto& iv : results) {
      shapely::Score& score = iv.second.front();
      raw_out << std::setw(4) << score.GetRotationType() << ": "
              << std::setw(10) << score.GetRawIndex() << " ("
              << std::setw(3) << score.GetXIndex() << ", "
              << std::setw(3) << score.GetYIndex() << ", "
              << std::setw(3) << score.GetZIndex() << ")  "
              << std::setprecision(6) << std::fixed << score.GetScore()
              << " : " << test_files[iv.first].stem() << std::endl;
      mat_out << std::setw(10) << score.GetScore();
      ss.str("");
      ss << "/Users/ivanwelsh/GitHub/shapely-drugs/data/ESR_DRUGS/TestOverlays/";
      ss << p1.stem().native() << ".test." << test_files[iv.first].stem().native() << "." << std::setw(2) << std::setfill('0') << score.GetRotationType() << ".xyz";
      std::ofstream test_out(ss.str());
      runner.PrintScoredTestXYZ(test_out, iv.first, score);
      test_out.close();
    }
    
    mat_out << "  " << p1.stem() << std::endl;
    
  }
  
  
  return 0;
}

int main(int argc, const char * argv[]) {
  return main_fftw_allbyall(argc, argv);
//  return main_esgrid_test(argc, argv);
// return file_check_main(argc, argv); 
//  return main_rdkittest(argc, argv);
  
#ifdef DEBUG
//  omp_set_num_threads(1);
#endif
  
//  return print_coordinates(argc, argv);
  
//  prep_main(argc, argv);
//  return scoring_main(argc, argv);
}



