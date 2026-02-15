#pragma once
#include <fstream>
#include <string>

class DataLogger {
    std::ofstream file;

public:
    DataLogger(const std::string& filename) {
        file.open(filename);
        file << "# step E dE |P| d|P| |L| d|L|\n";
    }

    void log(
        int step,
        double E, double dE,
        double P, double dP,
        double L, double dL
    ) {
        file << step << ","
             << E << "," << dE << ","
             << P << "," << dP << ","
             << L << "," << dL << "\n";
    }
};

// DataLogger es una clase simple para registrar datos de simulación en un archivo CSV.
// El constructor abre un archivo y escribe una línea de encabezado.