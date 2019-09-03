#include "systemcmd.h"

namespace sc
{
    SystemCmd::SystemCmd(const Parameter& parameter) : System(parameter)
    {

    }

    void SystemCmd::progressCallback(int current_iteration)
    {
        m_progressbar.progress(current_iteration);
    }

    bool SystemCmd::tearDownRun()
    {
        m_progressbar.reset();
        return true;
    }

    bool SystemCmd::initRun()
    {
        m_progressbar.setTotalTicks(this->simTime());
        m_progressbar.setEnabled(this->parameter().logParameter().showProgressBar());
        return true;
    }
}
