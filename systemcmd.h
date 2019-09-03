#ifndef SYSTEMCMD_H
#define SYSTEMCMD_H

#include "system.h"
#include "consoleprogressbar.h"

namespace sc
{
    class SystemCmd : public System
    {
      public:
        SystemCmd(const Parameter& parameter);

      protected:
        void progressCallback(int current_iteration) override final;
        bool tearDownRun() override final;
        bool initRun() override final;

      private:
        ConsoleProgressBar m_progressbar;
    };
}
#endif  // SYSTEMCMD_H
