"""Python wrapper for libswash_engine (ctypes).

    from swash_engine import Swash
    sw = Swash("/path/to/build/lib/libswash_engine.so")
    sw.init("/path/to/run/dir")          # directory holding INPUT and its data files
    while sw.step(10) > 0:
        eta = sw.field("WATL")            # (ny, nx) float32, -99 at inactive points
        sw.set_pressure(p)                # (ny, nx) float32, Pa
    sw.finish()

See engine/swash_engine.h for the semantics of each call.
"""
import ctypes
import os

import numpy as np

EXC = -99.0


class SwashError(RuntimeError):
    pass


class Swash:
    def __init__(self, libpath=None):
        libpath = libpath or os.environ.get("SWASH_ENGINE_LIB", "libswash_engine.so")
        self._lib = ctypes.CDLL(libpath)
        L = self._lib
        L.swash_init.argtypes = [ctypes.c_char_p]
        L.swash_init.restype = ctypes.c_int
        L.swash_step.argtypes = [ctypes.c_int]
        L.swash_step.restype = ctypes.c_int
        L.swash_status.restype = ctypes.c_int
        L.swash_finish.restype = ctypes.c_int
        L.swash_time.restype = ctypes.c_double
        L.swash_dt.restype = ctypes.c_double
        L.swash_end_time.restype = ctypes.c_double
        L.swash_grid_info.argtypes = [ctypes.POINTER(ctypes.c_int)] * 4 + [ctypes.POINTER(ctypes.c_double)] * 4
        L.swash_get_field.argtypes = [ctypes.c_char_p, ctypes.POINTER(ctypes.c_float), ctypes.c_int]
        L.swash_get_field.restype = ctypes.c_int
        for fn in ("swash_set_pressure_field", "swash_set_mass_source", "swash_set_bathymetry"):
            getattr(L, fn).argtypes = [ctypes.POINTER(ctypes.c_float)]
            getattr(L, fn).restype = ctypes.c_int
        L.swash_clear_pressure_field.restype = ctypes.c_int
        L.swash_save_state.argtypes = [ctypes.c_int]
        L.swash_save_state.restype = ctypes.c_int
        L.swash_restore_state.argtypes = [ctypes.c_int]
        L.swash_restore_state.restype = ctypes.c_int
        L.swash_acc_enabled.restype = ctypes.c_int
        self._cb_type = ctypes.CFUNCTYPE(None, ctypes.c_double)
        L.swash_set_step_callback.argtypes = [self._cb_type]
        self._cb = None
        self.nx = self.ny = self.nk = self.npts = 0
        self.x0 = self.y0 = self.dx = self.dy = 0.0

    # ---- lifecycle
    def init(self, workdir):
        rc = self._lib.swash_init(os.fspath(workdir).encode())
        if rc < 0:
            raise SwashError(f"swash_init failed with code {rc} (check PRINT / Errfile in {workdir})")
        nx, ny, nk, npts = (ctypes.c_int() for _ in range(4))
        x0, y0, dx, dy = (ctypes.c_double() for _ in range(4))
        self._lib.swash_grid_info(nx, ny, nk, npts, x0, y0, dx, dy)
        self.nx, self.ny, self.nk, self.npts = nx.value, ny.value, nk.value, npts.value
        self.x0, self.y0, self.dx, self.dy = x0.value, y0.value, dx.value, dy.value
        return rc

    def step(self, n=1):
        rc = self._lib.swash_step(int(n))
        if rc < 0:
            raise SwashError("swash_step failed (see PRINT / Errfile)")
        return rc

    def finish(self):
        return self._lib.swash_finish()

    def set_step_callback(self, fn):
        """fn(t: float) is called after every completed step."""
        self._cb = self._cb_type(fn) if fn is not None else self._cb_type(0)
        self._lib.swash_set_step_callback(self._cb)

    # ---- queries
    @property
    def status(self):
        return self._lib.swash_status()

    @property
    def time(self):
        return self._lib.swash_time()

    @property
    def dt(self):
        return self._lib.swash_dt()

    @property
    def end_time(self):
        return self._lib.swash_end_time()

    @property
    def gpu(self):
        return bool(self._lib.swash_acc_enabled())

    def field(self, name, k=1, masked=False):
        out = np.empty((self.ny, self.nx), dtype=np.float32)
        rc = self._lib.swash_get_field(name.encode(), out.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), int(k))
        if rc != 0:
            raise SwashError(f"swash_get_field({name!r}, k={k}) returned {rc}")
        if masked:
            return np.ma.masked_values(out, EXC)
        return out

    # ---- forcing
    def _dense(self, arr):
        a = np.ascontiguousarray(arr, dtype=np.float32)
        if a.shape != (self.ny, self.nx):
            raise ValueError(f"expected shape {(self.ny, self.nx)}, got {a.shape}")
        return a

    def set_pressure(self, pa):
        a = self._dense(pa)
        rc = self._lib.swash_set_pressure_field(a.ctypes.data_as(ctypes.POINTER(ctypes.c_float)))
        if rc != 0:
            raise SwashError(f"swash_set_pressure_field returned {rc}")

    def clear_pressure(self):
        self._lib.swash_clear_pressure_field()

    def set_mass_source(self, srcm):
        a = self._dense(srcm)
        rc = self._lib.swash_set_mass_source(a.ctypes.data_as(ctypes.POINTER(ctypes.c_float)))
        if rc != 0:
            raise SwashError(f"swash_set_mass_source returned {rc}")

    def set_bathymetry(self, depth):
        a = self._dense(depth)
        rc = self._lib.swash_set_bathymetry(a.ctypes.data_as(ctypes.POINTER(ctypes.c_float)))
        if rc != 0:
            raise SwashError(f"swash_set_bathymetry returned {rc}")

    # ---- checkpoints
    def save(self, slot=1):
        rc = self._lib.swash_save_state(int(slot))
        if rc != 0:
            raise SwashError(f"swash_save_state returned {rc}")

    def restore(self, slot=1):
        rc = self._lib.swash_restore_state(int(slot))
        if rc != 0:
            raise SwashError(f"swash_restore_state returned {rc}")


def hull_pressure(sw, x, y, heading, length, beam, draft, rho=1025.0, g=9.81, shape="box"):
    """Pressure footprint [Pa] of a vessel on the dense grid: a translated and
    rotated rectangle (or a parabolic 'hull' profile) of rho*g*draft.  This is
    the moving-pressure-patch representation of a ship (Kelvin wake technique)."""
    X = sw.field("XP"); Y = sw.field("YP")
    c, s = np.cos(heading), np.sin(heading)
    xl = (X - x) * c + (Y - y) * s
    yl = -(X - x) * s + (Y - y) * c
    inside = (np.abs(xl) <= length / 2) & (np.abs(yl) <= beam / 2)
    if shape == "box":
        prof = np.ones_like(X)
    else:  # smooth hull: parabolic in both directions
        prof = (1 - (2 * xl / length) ** 2) * (1 - (2 * yl / beam) ** 2)
    p = np.where(inside, rho * g * draft * prof, 0.0)
    return p.astype(np.float32)
