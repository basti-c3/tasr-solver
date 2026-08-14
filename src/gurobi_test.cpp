#include <gurobi_c.h>

#include <iostream>

int main() {
    GRBenv* env = nullptr;

    int error = GRBloadenv(&env, "gurobi-test.log");

    if (error) {
        std::cerr << "GRBloadenv failed with error "
            << error << '\n';

        if (env) {
            std::cerr << GRBgeterrormsg(env) << '\n';
        }

        return 1;
    }

    std::cout << "Gurobi environment successfully created.\n";

    GRBfreeenv(env);
    return 0;
}