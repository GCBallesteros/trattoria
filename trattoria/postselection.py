from typing import List, Tuple

import numpy as np

import trattoria


def construct_postselect_vector(
    timetrace_result: trattoria.TimeTraceResult,
    threshold: float,
    above: bool = True,
) -> List[Tuple[int, int]]:
    """Constructs a postselection vector based on thresholding the intensity trace.

    Ranges of records were the intensity is above/below some threshold will be excluded.

    Arguments
    ----------
    timetrace_result
        Intensity trace the postselection ranges are going to be based on.
    threshold
        Number of photons per bin set us threshold.
    above
        If True, return record ranges where intensity > threshold. If False, return
        record ranges where intensity < threshold.

    Returns
    -------
    List of 2-tuples (start_index, stop_index) of array indices in timetrace_y for
    each post-selected range of points.
    """
    intensity = timetrace_result.tt
    recnum = timetrace_result.recnum

    if above:
        if intensity[0] > threshold:
            add_first_time = True
        else:
            add_first_time = False
    else:
        if intensity[0] < threshold:
            add_first_time = True
        else:
            add_first_time = False

    # Find the transitions between on and off states by looking at the derivative.
    if above:
        select_ranges = np.diff(intensity > threshold, prepend=intensity[0]).nonzero()[
            0
        ]
    else:
        select_ranges = np.diff(intensity < threshold, prepend=intensity[0]).nonzero()[
            0
        ]

    if add_first_time:
        select_ranges = np.insert(select_ranges, 0, 0).astype(int)
    else:
        select_ranges = select_ranges[1:]

    if len(select_ranges) % 2 != 0:
        select_ranges = np.append(select_ranges, len(intensity) - 1).astype(int)

    select_ranges = select_ranges.reshape((-1, 2))

    # The select ranges are inclusive of the first index and exclusive of the second
    # one. Since the recnum include the last record on the bin when selecting the start
    # we must look at the previous recnum. On the second index because select_ranges is
    # exlusive we need to subtract one because we want the previous recnum.
    # In both cases we need to subtract one.
    select_ranges -= 1

    if select_ranges[0, 0] == -1:
        # The -1 index is conceptually the start of the measurement. i.e. recnum 0
        start_on_zero_record = True
        select_ranges[0, 0] = 0
    else:
        start_on_zero_record = False

    recnum_post_select_ranges = [
        (recnum[post_select_range[0]], recnum[post_select_range[1]])
        for post_select_range in select_ranges
    ]

    if start_on_zero_record:
        recnum_post_select_ranges[0] = (0, recnum_post_select_ranges[0][1])

    # There is a potential bug in that the last record of bins may be selected when
    # it is below treshold and it should not. It's not a big deal in the big scheme
    # of things but something to look out for.

    return recnum_post_select_ranges
