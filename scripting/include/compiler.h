#ifndef ATG_ENGINE_SIM_COMPILER_H
#define ATG_ENGINE_SIM_COMPILER_H

#include "language_rules.h"

#include "engine_sim.h"
#include "piranha.h"

#include <string>
#include <vector>

namespace es_script {

    class Compiler {
    public:
        struct Output {
            Engine *engine = nullptr;
            Vehicle *vehicle = nullptr;
            Transmission *transmission = nullptr;
            Simulator::Parameters simulatorParameters;
            ApplicationSettings applicationSettings;

            std::vector<Function *> functions;
        };

    private:
        static Output *s_output;

    public:
        Compiler();
        ~Compiler();

        static Output *output();

        void initialize();
        // nfs-mac: resolve imports without the working directory (an absolute directory searched after the
        // importing file's own), and write the error log to a chosen file; the host process need not chdir.
        void addSearchPath(const std::string &path);
        void setErrorLogPath(const std::string &path) { m_errorLogPath = path; }
        bool compile(const piranha::IrPath &path);
        Output execute();
        void destroy();

    private:
        void printError(const piranha::CompilationError *err, std::ofstream &file) const;

    private:
        LanguageRules m_rules;
        piranha::Compiler *m_compiler;
        piranha::NodeProgram m_program;
        std::string m_errorLogPath = "error_log.log";
    };

} /* namespace es_script */

#endif /* ATG_ENGINE_SIM_COMPILER_H */