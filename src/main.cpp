/*
 * Metodología de la Programación: Kmer4
 * Curso 2023/2024
 */

/**
 * @file main.cpp
 * @author Silvia Acid Carrillo <acid@decsai.ugr.es>
 * @author Andrés Cano Utrera <acu@decsai.ugr.es>
 * @author Luis Castillo Vidal <L.Castillo@decsai.ugr.es>
 * @author Javier Martínez Baena <jbaena@ugr.es>
 * @author estudiante1: Hugo Pereira Brito <hugopbrito@correo.ugr.es>
 * @author estudiante2: Álvaro Berenguer Cobo <alvarobc@correo.ugr.es>
 *
 * Created on 17 November 2023, 12:45
 */

#include "../include/Profile.h"
#include <cstring>
#include <iostream>

using namespace std;

/**
 * Shows help about the use of this program in the given output stream
 * @param outputStream The output stream where the help will be shown (for
 * example, cout, cerr, etc)
 */
void showEnglishHelp(ostream &outputStream) {
  outputStream << "ERROR in Kmer4 parameters" << endl;
  outputStream << "Run with the following parameters:" << endl;
  outputStream
      << "kmer4 [-t min|max] <file1.prf> <file2.prf> [ ... <filen.prf>]"
      << endl;
  outputStream << endl;
  outputStream << "Parameters:" << endl;
  outputStream << "-t min | -t max: search for minimun distances or maximum "
                  "distances (-t min by default)"
               << endl;
  outputStream << "<file1.prf>: source profile file for computing distances"
               << endl;
  outputStream << "<file2.prf> [ ... <filen.prf>]: target profile files for "
                  "computing distances"
               << endl;
  outputStream << endl;
  outputStream << "This program computes the distance from profile <file1.prf> "
                  "to the rest"
               << endl;
  outputStream << endl;
}

/**
 * This program reads an undefined number of Profile objects from the set of
 * files passed as parameters to main(). All the Profiles object, except the
 * first one, must be stored in a dynamic array of Profile objects. Then,
 * for each Profile in the dynamic array, this program prints to the
 * standard output the name of the file of that Profile and the distance from
 * the first Profile to the current Profile.
 * Finally, the program should print in the standard output, the name of
 * the file with the Profile with the minimum|maximum  distance to the Profile
 * of the first file and its profile identifier.
 *
 * At least, two Profile files are required to run this program.
 *
 * This program assumes that the profile files are already normalized and
 * sorted by frequency. This is not checked in this program. Unexpected results
 * will be obtained if those conditions are not met.
 *
 * Running sintax:
 * > kmer4 [-t min|max] <file1.prf> <file2.prf> [  ... <filen.prf>]
 *
 * Running example:
 * > kmer4 ../Genomes/human1.prf ../Genomes/worm1.prf ../Genomes/mouse1.prf
Distance to ../Genomes/worm1.prf: 0.330618
Distance to ../Genomes/mouse1.prf: 0.224901
Nearest profile file: ../Genomes/mouse1.prf
Identifier of the nearest profile: mus musculus
 *
 * Running example:
 * > kmer4 -t max ../Genomes/human1.prf ../Genomes/worm1.prf
../Genomes/mouse1.prf Distance to ../Genomes/worm1.prf: 0.330618 Distance to
../Genomes/mouse1.prf: 0.224901 Farthest profile file: ../Genomes/worm1.prf
Identifier of the farthest profile: worm
 */
int main(int argc, char *argv[]) {
  int util = argc - 2;

  // Process the main() arguments
  if (argc < 3) {
    showEnglishHelp(cerr);
    return 1;
  }

  if (*argv[1] == '-') {
    if (strcmp(argv[1], "-t") != 0) {
      showEnglishHelp(cerr);
      return 1;
    }
    if ((strcmp(argv[2], "min") != 0) && (strcmp(argv[2], "max")) != 0) {
      showEnglishHelp(cerr);
      return 1;
    }
    util = argc - 4;
  }

  // Allocate a dynamic array of Profiles
  Profile *profiles = new Profile[util];
  double *distances = new double[util];

  // Read the profiles
  Profile profile;
  profile.load(argv[argc - util - 1]);
  for (int i = 0; i < util; i++) {
    profiles[i].load(argv[i + argc - util]);
  }

  // Calculate and print the distance from the first Profile to the rest
  for (int i = 0; i < util; i++) {
    distances[i] = profile.getDistance(profiles[i]);
    cout << "Distance to " << argv[i + argc - util] << ": " << distances[i]
         << endl;
  }
  // Print name of the file and identifier that takes min|max distance to the
  // first one
  int profile_index = 0;
  int file_index = argc - util;
  if (strcmp(argv[1], "-t") == 0 && strcmp(argv[2], "max") == 0) {
    double max = distances[0];

    for (int i = 0; i < util; i++) {
      if (distances[i] > max) {
        max = distances[i];
        file_index = i + argc - util;
        profile_index = i;
      }
    }
    cout << "Farthest profile file: " << argv[file_index] << endl;
    cout << "Identifier of the farthest profile: "
         << profiles[profile_index].getProfileId() << endl;
  } else {
    double min = distances[0];

    for (int i = 0; i < util; i++) {
      if (distances[i] < min) {
        min = distances[i];
        file_index = i + argc - util;
        profile_index = i;
      }
    }
    cout << "Nearest profile file: " << argv[file_index] << endl;
    cout << "Identifier of the nearest profile: "
         << profiles[profile_index].getProfileId() << endl;
  }

  cout << profile.getDistance(profile);

  // Deallocate the dynamic array of Profile
  delete[] profiles;
  delete[] distances;
  return 0;
}
