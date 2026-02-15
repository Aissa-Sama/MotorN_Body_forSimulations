// io/regime_logger.h
#pragma once

#include <fstream>
#include <string>

// RegimeEvent ahora definido AQUÍ (ya no hay regime_event.h)
struct RegimeEvent {
    int step;
    int body_i;
    int body_j;
    std::string type;  // "ENTER_KS" | "EXIT_KS"
};

class RegimeLogger {
public:
    explicit RegimeLogger(const std::string& filename);
    ~RegimeLogger();

    void log(const RegimeEvent& event);

private:
    std::ofstream file;
};
