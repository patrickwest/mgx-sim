from io import StringIO

def window_contig(contig, max_window_size=5000, min_window_size=None):
    """
    Generate windowed contigs of length [max_window_size] from a reference contig. Returns a list of strings
    Default max_window_size: 5000
    Default min_window_size: 4999
    """
    
    # min window size default to max_window_size-1 (important for the end of contigs, can lower if you want to catch the shorter left over seqeunces)
    if min_window_size is None:
        min_window_size = int(max_window_size-1)
        
    # split contig into windows
    split_seqs = []
    s = StringIO(str(contig))
    while True:
        chunk = s.read(max_window_size)
        if len(chunk) > min_window_size:
            split_seqs.append(chunk)
        else:
            break
            
    return split_seqs
