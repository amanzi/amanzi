! closing missinging DSO objects 
      subroutine probinit(init,name,namlen,problo,probhi)
      integer init,namlen
      integer name(namlen)
      integer untin, i
      real*8  problo(2), probhi(2)
      end

      subroutine bl_proffortfuncstart(str)
      character*(*) str
      end

      subroutine bl_proffortfuncstop(str)
      character*(*) str
      end



