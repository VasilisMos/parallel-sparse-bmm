function C = bmmfiltered(A,B,F)
    C = ( F.*(A*B) ) > 0;
  end