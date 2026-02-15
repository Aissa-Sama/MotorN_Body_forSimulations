// io/regime_logger.cpp
#include "regime_logger.h"

RegimeLogger::RegimeLogger(const std::string& filename)
    : file(filename)
{
    file << "step,body_i,body_j,event\n";
}

RegimeLogger::~RegimeLogger() {
    if (file.is_open()) {
        file.close();
    }
}

void RegimeLogger::log(const RegimeEvent& event) {
    file << event.step << ","
         << event.body_i << ","
         << event.body_j << ","
         << event.type << "\n";
}