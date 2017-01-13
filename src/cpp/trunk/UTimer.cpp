#include "UTimer.h"
#include <sstream>
#include <cstdio>


UTimer::UTimer()
{
    Reset();
}

UTimer::~UTimer()
{

}

TimeT
UTimer::GetTicks()
{
    return time(NULL) - mStart;
}

TimeT
UTimer::TicksPerSec()
{
    return 1;
}


TimeT
UTimer::ClockMax()
{
    time_t tv;
    int bits = sizeof(tv) << 3;
    tv = (clock_t) -1 ^ (1 << (bits - 1));

    return  TimeT(tv);
}

void
UTimer::Reset()
{
    mStart = time(NULL);
}

TimeT
UTimer::DurationInSeconds(TimeT inTicks )
{
    return inTicks/TicksPerSec() % 60;
}

TimeT
UTimer::DurationInMinutes(TimeT inTicks)
{
    return inTicks/(TicksPerSec() * 60) % 60;
}

TimeT
UTimer::DurationInHours(TimeT inTicks)
{
    return inTicks/(TicksPerSec() * 3600) % 24;
}

TimeT
UTimer::DurationInDays(TimeT inTicks)
{
    return inTicks/(TicksPerSec() * 3600 * 24);
}

void
UTimer::Duration(TimeT inTicks, TimeT & outDays, TimeT & outHours, TimeT & outMinutes, TimeT & outSeconds)
{
    outSeconds = DurationInSeconds(inTicks);
    outMinutes = DurationInMinutes(inTicks);
    outHours = DurationInHours(inTicks);
    outDays = DurationInDays(inTicks);
}

string
UTimer::DurationInText(TimeT inTicks)
{
    char s[255];
    ostringstream ss;
    
    TimeT days, hours, minutes, seconds;

    Duration(inTicks, days, hours, minutes, seconds);

    if (days == 1)
        ss << days << " day and ";
    else if (days > 1)
        ss << days << " days and ";
        
    ss.fill('0');
    ss.width(2);
    ss <<   hours << ":";
    ss.width(2);
    ss << minutes << ":";
    ss.width(2);
    ss << seconds;

    return ss.str();
}

