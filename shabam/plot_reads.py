
from __future__ import division

from shabam.colors import COLORS

def plot_read(context, bases, quals=None, x_offset=0, y_offset=0, width=None,
        is_reverse=False, by_strand=False, merge_delta=0.05):
    ''' plots the bases in a read to a cairocffi.Context
    
    Args:
        context: cairocffi.Context as a plotting device
        bases: list of bases (per parse_read, so indels are odd)
        quals: list of quality scores for each base
        x_offset: x position to start plotting the read at
        y_offset: y position to plot the read at
        width: with of image in pixels
        is_reverse: whether the read is for the reverse strand
        by_strand: boolean for whether we want to shade reads by strand
        merge_delta: difference allowed between neighboring qualities before not merging
    '''
    
    if quals is None:
        quals = [100] * len(bases)
    
    if width is None:
        width = len(bases) * 10
    
    i = 0
    while i < len(bases):
        base = bases[i]
        qual = quals[i]
        x_pos = (x_offset + i) * 10
        
        # collapse adjacent bases if they have the same base call and quality
        delta = 1
        next_base = bases[i + delta] if (i + delta) < len(bases) else None
        next_qual = quals[i + delta] if (i + delta) < len(bases) else None
        
        # account for deleted bases with qual of '-'
        qual = 100 if isinstance(qual, str) else qual
        next_qual = qual if isinstance(next_qual, str) else next_qual
        
        while (base == next_base) and (abs(1 - relative_qual(qual, next_qual)) <  merge_delta):
            next_base = bases[i + delta] if (i + delta) < len(bases) else None
            next_qual = quals[i + delta] if (i + delta) < len(bases) else None
            next_qual = qual if isinstance(next_qual, str) else next_qual
            delta += 1
        
        i += delta
        
        if x_pos < 0 and x_pos + delta * 10 > 0:
            delta -= abs(x_pos) / 10
            x_pos = 0
        
        if x_pos < 0 or x_pos > width - 1:
            # don't plot bases outside the required window. This is necessary
            # when plotting SVGs, otherwise the SVG includes the outside bases.
            continue
        
        if len(base) > 1:
            plot_insertion(context, base, x_pos, y_offset)
            base = 'M'
            qual = qual[0]
        
        if base == 'M' and by_strand:
            strand = {True: 'r', False: 'f'}[is_reverse]
            base = f'M_{strand}'
        
        context.rectangle(x_pos, y_offset, 10 * delta, 10)
        context.set_source_rgba(*(COLORS[base.upper()] + [to_alpha(qual)]))
        context.fill()

def to_alpha(qual, threshold=35):
    ''' convert base quality to an alpha transparency float
    '''
    try:
        return min(threshold, qual)/threshold
    except TypeError:
        return 1.0

def relative_qual(qual: int, next_qual: int):
    ''' get the relative quality between two quality scores
    
    Has to handle cases where one or both scores are zero
    
    Args:
        qual: 
        next_qual:
    
    Returns:
        ratio of scores to each other
    '''
    if qual == 0:
        return 0 if (next_qual == 0 and qual == 0) else 1
    else:
        return next_qual / qual

def plot_insertion(context, bases, x_pos, y_offset):
    ''' plot inserted bases at the insertion site
    
    Args:
        context: cairocffi.Context as a plotting device
        bases: string of inserted bases
        x_pos: position of insertion (in pixels)
        y_offset: y position to plot the read at
    '''
    
    # select a font, and figure out the text sizes, so we can align text
    context.select_font_face('Arial')
    context.set_font_size(7)
    _, _, width, _, _, _ = context.text_extents(bases)
    
    context.move_to(x_pos + 10 - width/2, y_offset - 3)
    context.set_source_rgb(0, 0, 0)
    context.show_text(bases)
    
    # plot an arrow to indicate the insertion point
    context.set_line_width(1)
    context.move_to(x_pos + 10 - 1.5, y_offset - 1.5)
    context.line_to(x_pos + 10, y_offset)
    context.line_to(x_pos + 10 + 1.5, y_offset - 1.5)
    context.close_path()
    context.stroke_preserve()
    context.set_source_rgb(0, 0, 0)
    context.fill()
