// {"means":[2.0,2.0],"weights":[0.5,0.5],"stdevs":[0.0,0.0],"modal_n":2,"n_reads":[4,2],"seqs":[["GTAGTAG","single"],["GTAGTAG","single"]],"start_anchor_seqs":[["CTGCA","single"],["CTGCA","single"]]},"ps":1,"read_peaks_called":true,"time":0.04344633303117007}

export type Peaks = {
  means: number[];
  weights: number[];
  stdevs: number[];
  modal_n: number;
  n_reads: number[];
  seqs?: ([string, "single" | "best_rep" | "poa"])[],
  ps?: number;
  read_peaks_called?: boolean;
  time?: number;
};

export type Read = {

};
