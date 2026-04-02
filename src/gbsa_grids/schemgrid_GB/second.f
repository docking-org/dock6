c---------------------------------------------------------------------
c
c  --contains system-specific code
c  --current operating systems are VMS and UNIX
c  One of these options may be able to substitute for "second.dummy.f"
c  on your system.  Edit comp.com accordingly to try it out.
c
c  ECM  5/91, code from ARLeach.
c---------------------------------------------------------------------
      subroutine second(time)
c  --finds the cpu time used 
c
      dimension tarray(2)
c
c VMS specific code:
c       parameter jpi$c_listend = '00000000'x
c       parameter jpi$_cputim = '00000407'x
c       integer*4 i4blk(4)
c       integer*2 i2blk(8)
c       equivalence ( i4blk(1) , i2blk(1) )
c       data i2blk / 4, jpi$_cputim, 0, 0, 0, 0, jpi$c_listend, 0  /
c       i4blk(2) = %loc ( icpu )
c       istat = sys$getjpi ( , , ,  i2blk , , ,  )
c       if ( .not. istat ) rtime=0.0
c       time = float ( icpu ) / 100.
c
c UNIX specific code:
c (note:  on some machines, etime is double precision, so the following
c should be declared accordingly...)
        real time, etime, tarray
        time=etime(tarray)
c
      return
      end
c---------------------------------------------------------------------
