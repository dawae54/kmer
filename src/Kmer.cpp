/*
 * Metodología de la Programación: Kmer3
 * Curso 2023/2024
 */

/**
 * @file Kmer.cpp
 * @author Silvia Acid Carrillo <acid@decsai.ugr.es>
 * @author Andrés Cano Utrera <acu@decsai.ugr.es>
 * @author Luis Castillo Vidal <L.Castillo@decsai.ugr.es>
 * @author Javier Martínez Baena <jbaena@ugr.es>
 * @author estudiante1: Hugo Pereira Brito <hugopbrito@correo.ugr.es>
 * @author estudiante2: Álvaro Berenguer Cobo <alvarobc@correo.ugr.es>
 *
 * Created on 24 October 2023, 14:00
 */

#include "../include/Kmer.h"
#include <string>
using namespace std;

Kmer::Kmer(int k) : _text(k, MISSING_NUCLEOTIDE) {
    if (k <= 0) {
        throw std::invalid_argument("Number of nucleotides cannot be 0");
    }
}

Kmer::Kmer(const string &text) : _text(text) {
    if (text.size() < 1) {
        throw std::invalid_argument("The text cannot be empty");
    }
}

int Kmer::getK() const { return size(); }

int Kmer::size() const { return _text.size(); }

string Kmer::toString() const { return _text; }

const char &Kmer::at(int index) const { return _text.at(index); }

char &Kmer::at(int index) { return _text.at(index); }

void Kmer::toLower() { boost::algorithm::to_lower(this->_text); }

void Kmer::toUpper() { boost::algorithm::to_upper(this->_text); }

void Kmer::normalize(const std::string &validNucleotides) {
    // Pasamos a mayúsculas
    ToUpper(*this);

    // Reemplazamos por carácter nulo
    for (size_t i = 0; i < _text.size(); i++) {
        if (!IsValidNucleotide(_text[i], validNucleotides)) {
            _text.at(i) = '_';
        }
    }
}

Kmer Kmer::complementary(const string &nucleotides,
                         const string &complementaryNucleotides) const {
    if (nucleotides.size() != complementaryNucleotides.size())
        throw std::invalid_argument(
            "Error: \"nucleotides and \"complementaryNucleotides\" must have "
            "the same length");
    Kmer comp = *this;

    for (int i = 0; i < size(); i++) {
        if (IsValidNucleotide(at(i), nucleotides)) {
            comp.at(i) = complementaryNucleotides.at(nucleotides.find(at(i)));
        }
    }

    return comp;
}

int Kmer::find(char to_find) const { return _text.find(to_find); }

bool IsValidNucleotide(char nucleotide, const string &validNucleotides) {
    return validNucleotides.find(nucleotide) != string::npos;
}

void ToLower(Kmer &kmer) { kmer.toLower(); }

void ToUpper(Kmer &kmer) { kmer.toUpper(); }
