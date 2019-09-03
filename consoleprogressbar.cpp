#include "consoleprogressbar.h"
#include <algorithm>
#include <stdio.h>
#include <string>

namespace sc
{
    ConsoleProgressBar::ConsoleProgressBar()
        : m_enabled(true)
        , m_total_ticks(100)
        , m_width(72)
        , m_t_width(0)
        , m_last_perc(-1)
        , m_label("Progress ")
    {
        this->precalc();
    }

    void ConsoleProgressBar::progress(int tick)
    {
        if (!m_enabled)
        {
            return;
        }

        int percent = std::min(100, (tick * 100) / m_total_ticks);
        if (percent <= m_last_perc)
        {
            return;
        }

        // minus label len
        int pos = std::min(m_t_width, (tick * m_t_width) / m_total_ticks);

        printf("%s[", m_label.c_str());

        for (int i = 0; i < pos; i++)
        {
            printf("%c", '=');
        }

        printf("%*c", m_t_width - pos + 1, ']');
        printf(" %3d%%\r", percent);

        m_last_perc = percent;
    }

    void ConsoleProgressBar::reset()
    {
        m_last_perc = -1;
    }

    int ConsoleProgressBar::totalTicks() const
    {
        return m_total_ticks;
    }

    void ConsoleProgressBar::setTotalTicks(int total_ticks)
    {
        if (total_ticks >= 0)
        {
            m_total_ticks = total_ticks;
        }
    }

    int ConsoleProgressBar::width() const
    {
        return m_width;
    }

    void ConsoleProgressBar::setWidth(int width)
    {
        m_width = width;
        this->precalc();
    }

    std::string ConsoleProgressBar::label() const
    {
        return m_label;
    }

    void ConsoleProgressBar::setLabel(const std::string& label)
    {
        m_label = label;
        this->precalc();
    }

    void ConsoleProgressBar::precalc()
    {
        m_t_width = m_width - m_label.length();
    }

    bool ConsoleProgressBar::enabled() const
    {
        return m_enabled;
    }

    void ConsoleProgressBar::setEnabled(bool enabled)
    {
        m_enabled = enabled;
    }
}
