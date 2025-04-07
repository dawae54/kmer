/*
 * Metodología de la Programación: Kmer3
 * Curso 2023/2024
 */

/**
 * @file Profile.cpp
 * @author Silvia Acid Carrillo <acid@decsai.ugr.es>
 * @author Andrés Cano Utrera <acu@decsai.ugr.es>
 * @author Luis Castillo Vidal <L.Castillo@decsai.ugr.es>
 * @author Javier Martínez Baena <jbaena@ugr.es>
 * @author estudiante1: Hugo Pereira Brito <hugopbrito@correo.ugr.es>
 * @author estudiante2: Álvaro Berenguer Cobo <alvarobc@correo.ugr.es>
 *
 * Created on 29 January 2023, 11:00
 */

#include "../include/Profile.h"
#include <cmath>

using namespace std;

const string Profile::MAGIC_STRING_T = "MP-KMER-T-1.0";

Profile::Profile() {
    this->_profileId = "unknown";
    this->_vectorKmerFreq.reserve(INITIAL_CAPACITY);
}

Profile::Profile(int size) {
    this->_profileId = "unknown";
    if (size <= 0)
        throw std::out_of_range("Size must be greater than 0");

    this->_vectorKmerFreq.reserve(size);
}

Profile::Profile(const Profile &orig) {
    this->_profileId = orig.getProfileId();
    this->_vectorKmerFreq = orig._vectorKmerFreq;
}

Profile::~Profile() {
    // deallocate();
}

Profile Profile::operator=(const Profile &orig) {
    if (this->_vectorKmerFreq.data() != orig._vectorKmerFreq.data()) {
        this->_vectorKmerFreq = orig._vectorKmerFreq;
        this->_profileId = orig._profileId;
    }
    return *this;
}

string Profile::getProfileId() const { return this->_profileId; }

void Profile::setProfileId(const std::string &id) { this->_profileId = id; }

const KmerFreq &Profile::at(int index) const {
    return this->_vectorKmerFreq.at(index);
}

KmerFreq &Profile::at(int index) { return this->_vectorKmerFreq.at(index); }

int Profile::getSize() const { return this->_vectorKmerFreq.size(); }

int Profile::getCapacity() const { return this->_vectorKmerFreq.capacity(); }

double Profile::getDistance(const Profile &otherProfile) const {
    if (this->_vectorKmerFreq.empty() || otherProfile._vectorKmerFreq.empty())
        throw ::std::invalid_argument(
            "Every Profile must have at least one KmerFreq");
    double distance = 0;
    int rank;

    for (int i = 0; i < getSize(); i++) {
        rank = otherProfile.findKmer(at(i).getKmer());
        if (rank == -1)
            rank = otherProfile.getSize();
        distance = distance + abs(i - rank);
    }

    distance = distance / (getSize() * otherProfile.getSize());
    return distance;
}

int Profile::findKmer(const Kmer &kmer, int initialPos, int finalPos) const {
    bool found = false;
    int foundPos = -1;
    for (int i = initialPos; i <= finalPos && found == false; i++) {
        if (kmer.toString() == this->_vectorKmerFreq[i].getKmer().toString()) {
            foundPos = i;
            found = true;
        }
    }

    return foundPos;
}

int Profile::findKmer(const Kmer &kmer) const {
    //    int foundPos = -1;
    //    for (int i=0; i<this->getSize(); i++){
    //        if (kmer.toString() == _vectorKmerFreq[i].getKmer().toString())
    //            foundPos = i;
    //    }
    return findKmer(kmer, 0, getSize() - 1);
    //    return foundPos;
}

string Profile::toString() const {
    string output;
    output += this->_profileId + "\n";
    output += to_string(this->getSize()) + "\n";

    for (int i = 0; i < this->getSize(); i++) {
        output += this->_vectorKmerFreq[i].toString() + "\n";
    }

    return output;
}

void Profile::sort() {
    KmerFreq max;
    int max_pos;

    for (int i = 0; i < this->getSize(); i++) {
        max = _vectorKmerFreq[i];
        max_pos = i;
        for (int j = i + 1; j < this->getSize(); j++) {
            if (_vectorKmerFreq[j].getFrequency() > max.getFrequency()) {
                max = _vectorKmerFreq[j];
                max_pos = j;
            } else if (_vectorKmerFreq[j].getFrequency() ==
                           max.getFrequency() &&
                       _vectorKmerFreq[j].getKmer().toString() <
                           max.getKmer().toString()) {
                max = _vectorKmerFreq[j];
                max_pos = j;
            }
        }

        _vectorKmerFreq[max_pos] = _vectorKmerFreq[i];
        _vectorKmerFreq[i] = max;
    }
}

void Profile::save(const char *fileName) const {
    ofstream file;
    // Abrimos fichero
    file.open(fileName);
    if (!file.is_open())
        throw std::ios_base::failure("Unable to open file");
    file << MAGIC_STRING_T << "\n";
    file << toString();

    file.close();
}

void Profile::load(const char *fileName) {
    // Abrimos fichero
    ifstream file;
    file.open(fileName);
    if (!file.is_open())
        throw std::ios_base::failure("Unable to open file");

    string magic_string;
    getline(file, magic_string);
    if (magic_string != MAGIC_STRING_T)
        throw std::invalid_argument("Invalid magic string");

    getline(file, this->_profileId);

    int size;
    file >> size;
    // file.ignore();

    if (size <= 0)
        throw std::out_of_range("Size must be greater or equal to 0");

    // Borramos lo que había
    _vectorKmerFreq.clear();

    // Reservamos
    _vectorKmerFreq.reserve(size);

    // Guardamos los datos del fichero
    for (int i = 0; i < size; i++) {
        KmerFreq kmer_freq;
        string kmer;
        int frequency;

        file >> kmer;
        file >> frequency;

        kmer_freq.setKmer(Kmer(kmer));
        kmer_freq.setFrequency(frequency);
        this->append(kmer_freq);
    }

    file.close();
}

void Profile::append(const KmerFreq &kmerFreq) {
    int found = findKmer(kmerFreq.getKmer());
    if (found != -1) {
        this->_vectorKmerFreq[found].setFrequency(
            this->_vectorKmerFreq[found].getFrequency() + 1);
    }
}

void Profile::normalize(const std::string &validNucleotides) {
    int found;
    Kmer normalized;
    // Loop to traverse and normalize each one of the kmers in array
    // Normalize kmer i
    for (int i = 0; i < this->getSize(); i++) {
        normalized = _vectorKmerFreq[i].getKmer();
        normalized.normalize(validNucleotides);
        _vectorKmerFreq[i].setKmer(normalized);
    }

    for (int i = 1; i < this->getSize();) {
        found = findKmer(_vectorKmerFreq[i].getKmer(), 0, i - 1);
        if (found >= 0) {
            _vectorKmerFreq[found].setFrequency(
                _vectorKmerFreq[found].getFrequency() +
                _vectorKmerFreq[i].getFrequency());
            deletePos(i);
        } else
            i++;
    }
}

void Profile::deletePos(int pos) {
    this->_vectorKmerFreq.erase(this->_vectorKmerFreq.begin() + pos);
}

void Profile::zip(bool deleteMissing, int lowerBound) {
    for (int i = 0; i < this->getSize(); i++) {
        if ((_vectorKmerFreq[i].getFrequency() <= lowerBound) ||
            (deleteMissing &&
             _vectorKmerFreq[i].getKmer().find(
                 _vectorKmerFreq[i].getKmer().MISSING_NUCLEOTIDE) != -1)) {
            deletePos(i);
            i--;
        }
    }
}

void Profile::join(const Profile &profile) {
    for (int i = 0; i < profile.getSize(); i++) {
        this->append(profile.at(i));
    }
}

// void Profile::allocate(int capacity)
// {
//     if (capacity > 0)
//     {
//         _vectorKmerFreq = new KmerFreq[capacity];
//         _capacity = capacity;
//         _size = 0;
//     }
// }

// void Profile::reallocate(int capacity)
// {
//     if (capacity <= 0)
//     {
//         throw std::out_of_range("Capacity must be greater than 0");
//     }
//     if (capacity > _capacity)
//     {
//         Profile aux;
//         aux.allocate(capacity);
//         aux.copy(*this);
//         this->deallocate();
//         this->_vectorKmerFreq = aux._vectorKmerFreq;
//     }
//     else
//     {
//         deallocate();
//         allocate(capacity);
//     }
// }

// void Profile::deallocate()
// {
//     if (_capacity > 0)
//     {
//         delete[] _vectorKmerFreq;
//         _vectorKmerFreq = nullptr;
//         _capacity = 0;
//         _size = 0;
//     }
// }

// void Profile::copy(const Profile &to_copy)
// {
//     if (to_copy.getSize() <= this->_capacity)
//     {
//         for (int i = 0; i < to_copy.getSize(); i++)
//         {
//             this->_vectorKmerFreq[i] = to_copy.at(i);
//         }
//         this->_size = to_copy.getSize();
//     }
// }
