#ifndef SC_CONSOLEPROGRESSBAR_H
#define SC_CONSOLEPROGRESSBAR_H

#include <string>

namespace sc
{
    class ConsoleProgressBar
    {
      public:
        ConsoleProgressBar();
        void progress(int tick);
        void reset();

        int totalTicks() const;
        void setTotalTicks(int total_ticks);

        int width() const;
        void setWidth(int width);

        std::string label() const;
        void setLabel(const std::string &label);

        bool enabled() const;
        void setEnabled(bool enabled);

    private:
        void precalc();

    private:
        bool m_enabled;
        int m_total_ticks;
        int m_width;
        int m_t_width;
        int m_last_perc;
        std::string m_label;
    };
}

#endif  // CONSOLEPROGRESSBAR_H
