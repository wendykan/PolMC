#ifndef __UTIMER_H
#define __UTIMER_H

#include <ctime>
#include <string>

using namespace std;

typedef time_t TimeT;

class UTimer {
public:
    UTimer();
    virtual ~UTimer();

    TimeT GetTicks();
    TimeT TicksPerSec();
    TimeT ClockMax();
    void Reset();
    TimeT	DurationInSeconds(TimeT inTicks );
    TimeT	DurationInMinutes(TimeT inTicks );
    TimeT	DurationInHours(TimeT inTicks );
    TimeT	DurationInDays(TimeT inTicks );
    void	Duration(TimeT inTicks,
                   TimeT & outDays,
                   TimeT & outHours,
                   TimeT & outMinutes,
                   TimeT & outSeconds);

    string	DurationInText(TimeT inTicks);

protected:
        TimeT	mUnwrappedTimer;
    TimeT	mWrap;
    time_t			mStart;
    time_t			mPrevious;
};

#endif
