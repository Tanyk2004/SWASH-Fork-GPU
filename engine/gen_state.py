#!/usr/bin/env python3
"""Regenerates engine/swash_state.ftn90 from the allocatable arrays declared in
src/SwashFlowdata.ftn90 and src/SwashRigBoddata.ftn90.

Run from the repository root after changing either module:

    python3 engine/gen_state.py
"""
import pathlib
import re

root = pathlib.Path(__file__).resolve().parents[1]
pat = re.compile(r'^ *([a-z*0-9]+) *, *dimension\(([^)]*)\) *, *save *, *allocatable *:: *([A-Za-z0-9_]+)', re.M)
flow = [m for m in pat.findall((root / 'src/SwashFlowdata.ftn90').read_text()) if not m[0].startswith('type')]
rig = pat.findall((root / 'src/SwashRigBoddata.ftn90').read_text())


def comp(t, d, n):
    return f"       {t:<8s}, dimension({d}), allocatable :: {n}\n"


def save(t, d, n):
    r = d.count(':')
    b = ",".join(f"lbound({n},{i}):ubound({n},{i})" for i in range(1, r + 1))
    return (f"       if ( allocated({n}) ) then\n"
            f"          if ( allocated(sl%{n}) ) deallocate(sl%{n})\n"
            f"          allocate(sl%{n}({b}))\n"
            f"          sl%{n} = {n}\n"
            f"       endif\n")


def rest(t, d, n):
    return f"       if ( allocated(sl%{n}) .and. allocated({n}) ) {n} = sl%{n}\n"


decl = "".join(comp(*a) for a in flow) + "       ! rigid-body state\n" + "".join(comp(*a) for a in rig)
sv = "".join(save(*a) for a in flow) + "       !\n       ! rigid-body state\n       !\n" + "".join(save(*a) for a in rig)
rs = "".join(rest(*a) for a in flow) + "       !\n" + "".join(rest(*a) for a in rig)

out = root / 'engine/swash_state.ftn90'
text = out.read_text()
head = text[:text.index('    type state_slot')]
head = re.sub(r'\(\d+ \+ \d+ arrays\)', f'({len(flow)} + {len(rig)} arrays)', head)
body = f"""    type state_slot
       logical :: used = .false.
       real*8  :: timco
       real*8  :: dt
       real    :: cflmax
       integer :: it
       integer :: istep
       character(20) :: chtime
{decl}    end type state_slot
!
    type(state_slot), dimension(max_slots), save, target :: slots
!
contains
!
    subroutine state_save ( islot )
       implicit none
       integer, intent(in) :: islot
       type(state_slot), pointer :: sl
       !
       sl => slots(islot)
       sl%used   = .true.
       sl%timco  = timco
       sl%dt     = dt
       sl%cflmax = cflmax
       sl%it     = it
       sl%istep  = istep
       sl%chtime = chtime
       !
{sv}    end subroutine state_save
!
    subroutine state_restore ( islot )
       implicit none
       integer, intent(in) :: islot
       type(state_slot), pointer :: sl
       !
       sl => slots(islot)
       if ( .not.sl%used ) return
       timco  = sl%timco
       dt     = sl%dt
       cflmax = sl%cflmax
       it     = sl%it
       istep  = sl%istep
       chtime = sl%chtime
       !
{rs}    end subroutine state_restore
!
    logical function state_used ( islot )
       implicit none
       integer, intent(in) :: islot
       state_used = slots(islot)%used
    end function state_used
!
end module swash_state
"""
out.write_text(head + body)
print(f"wrote {out} with {len(flow)} + {len(rig)} arrays")
