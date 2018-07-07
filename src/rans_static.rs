type c_void = u8;
type c_char = u8;
type c_uchar = u8;
type c_schar = u8;
type c_ulong = u64; // sorry, it's inconsistent
type c_long = i64;
type c_ulonglong = u64;
type c_longlong = i64;
type c_ushort = u16;
type c_short = i16;
type c_uint = u32;
type c_int = i32;
type c_double = f64;
type c_float = f32;
extern crate libc;
extern "C" {
    pub type _IO_FILE_plus;
    #[no_mangle]
    static mut _IO_2_1_stdin_: _IO_FILE_plus;
    #[no_mangle]
    static mut _IO_2_1_stdout_: _IO_FILE_plus;
    #[no_mangle]
    static mut _IO_2_1_stderr_: _IO_FILE_plus;
    #[no_mangle]
    fn __uflow(_: *mut _IO_FILE) -> c_int;
    #[no_mangle]
    fn __overflow(_: *mut _IO_FILE, _: c_int) -> c_int;
    #[no_mangle]
    fn _IO_getc(__fp: *mut _IO_FILE) -> c_int;
    #[no_mangle]
    fn _IO_putc(__c: c_int, __fp: *mut _IO_FILE) -> c_int;
    #[no_mangle]
    static mut stdin: *mut _IO_FILE_0;
    #[no_mangle]
    static mut stdout: *mut _IO_FILE_0;
    #[no_mangle]
    static mut stderr: *mut _IO_FILE_0;
    #[no_mangle]
    fn fprintf(_: *mut FILE, _: *const c_char, ...) -> c_int;
    #[no_mangle]
    fn vfprintf(_: *mut FILE, _: *const c_char, _: *mut __va_list_tag)
     -> c_int;
    #[no_mangle]
    static mut sys_nerr: c_int;
    #[no_mangle]
    static sys_errlist: [*const c_char; 0];
    #[no_mangle]
    fn strtod(__nptr: *const c_char, __endptr: *mut *mut c_char)
     -> c_double;
    #[no_mangle]
    fn strtol(__nptr: *const c_char, __endptr: *mut *mut c_char,
              __base: c_int) -> c_long;
    #[no_mangle]
    fn strtoll(__nptr: *const c_char, __endptr: *mut *mut c_char,
               __base: c_int) -> c_longlong;
    #[no_mangle]
    fn malloc(_: c_ulong) -> *mut c_void;
    #[no_mangle]
    fn free(__ptr: *mut c_void) -> ();
    #[no_mangle]
    fn exit(_: c_int) -> !;
    #[no_mangle]
    static mut __environ: *mut *mut c_char;
    #[no_mangle]
    static mut optarg: *mut c_char;
    #[no_mangle]
    static mut optind: c_int;
    #[no_mangle]
    static mut opterr: c_int;
    #[no_mangle]
    static mut optopt: c_int;
    #[no_mangle]
    fn memmove(_: *mut c_void, _: *const c_void, _: c_ulong)
     -> *mut c_void;
    #[no_mangle]
    fn memset(_: *mut c_void, _: c_int, _: c_ulong)
     -> *mut c_void;
    #[no_mangle]
    fn memcmp(_: *const c_void, _: *const c_void,
              _: c_ulong) -> c_int;
    #[no_mangle]
    fn __rawmemchr(__s: *const c_void, __c: c_int)
     -> *mut c_void;
    #[no_mangle]
    fn gettimeofday(__tv: *mut timeval, __tz: __timezone_ptr_t)
     -> c_int;
}
pub type uint32_t = c_uint;
pub type __off_t = c_long;
pub type _IO_FILE = _IO_FILE_0;
pub type __off64_t = c_long;
#[derive ( Copy , Clone )]
#[repr ( C )]
pub struct _IO_FILE_0 {
    pub _flags: c_int,
    pub _IO_read_ptr: *mut c_char,
    pub _IO_read_end: *mut c_char,
    pub _IO_read_base: *mut c_char,
    pub _IO_write_base: *mut c_char,
    pub _IO_write_ptr: *mut c_char,
    pub _IO_write_end: *mut c_char,
    pub _IO_buf_base: *mut c_char,
    pub _IO_buf_end: *mut c_char,
    pub _IO_save_base: *mut c_char,
    pub _IO_backup_base: *mut c_char,
    pub _IO_save_end: *mut c_char,
    pub _markers: *mut _IO_marker,
    pub _chain: *mut _IO_FILE_0,
    pub _fileno: c_int,
    pub _flags2: c_int,
    pub _old_offset: __off_t,
    pub _cur_column: c_ushort,
    pub _vtable_offset: c_schar,
    pub _shortbuf: [c_char; 1],
    pub _lock: *mut c_void,
    pub _offset: __off64_t,
    pub __pad1: *mut c_void,
    pub __pad2: *mut c_void,
    pub __pad3: *mut c_void,
    pub __pad4: *mut c_void,
    pub __pad5: size_t,
    pub _mode: c_int,
    pub _unused2: [c_char; 20],
}
pub type Rans64State = uint64_t;
#[derive ( Copy , Clone )]
#[repr ( C )]
pub struct __va_list_tag {
    pub gp_offset: c_uint,
    pub fp_offset: c_uint,
    pub overflow_arg_area: *mut c_void,
    pub reg_save_area: *mut c_void,
}
#[derive ( Copy , Clone )]
#[repr ( C )]
pub struct ari_decoder {
    pub R: [c_uchar; 4096],
}
pub type _IO_lock_t = ();
pub type __compar_fn_t =
    Option<unsafe extern "C" fn(_: *const c_void,
                                _: *const c_void) -> c_int>;
#[derive ( Copy , Clone )]
#[repr ( C )]
pub struct timeval {
    pub tv_sec: __time_t,
    pub tv_usec: __suseconds_t,
}
#[derive ( Copy , Clone )]
#[repr ( C )]
pub struct RansEncSymbol {
    pub rcp_freq: uint64_t,
    pub freq: uint32_t,
    pub bias: uint32_t,
    pub cmpl_freq: uint32_t,
    pub rcp_shift: uint32_t,
}
#[derive ( Copy , Clone )]
#[repr ( C )]
pub struct _IO_marker {
    pub _next: *mut _IO_marker,
    pub _sbuf: *mut _IO_FILE_0,
    pub _pos: c_int,
}
pub type Rans64DecSymbol = RansDecSymbol;
#[derive ( Copy , Clone )]
#[repr ( C )]
pub struct RansDecSymbol {
    pub start: uint32_t,
    pub freq: uint32_t,
}
pub type FILE = _IO_FILE_0;
pub type size_t = c_ulong;
pub type __suseconds_t = c_long;
pub type __timezone_ptr_t = *mut timezone;
pub type RansState = Rans64State;
pub type __time_t = c_long;
#[derive ( Copy , Clone )]
#[repr ( C )]
pub struct timezone {
    pub tz_minuteswest: c_int,
    pub tz_dsttime: c_int,
}
pub type Rans64EncSymbol = RansEncSymbol;
pub type uint64_t = c_ulong;
pub type uint8_t = c_uchar;
#[derive ( Copy , Clone )]
#[repr ( C )]
pub struct blocks {
    pub blk: *mut c_uchar,
    pub sz: uint32_t,
}
unsafe extern "C" fn vprintf(mut __fmt: *const c_char,
                             mut __arg: *mut __va_list_tag) -> c_int {
    return vfprintf(stdout, __fmt, __arg);
}
unsafe extern "C" fn getchar() -> c_int { return _IO_getc(stdin); }
unsafe extern "C" fn getc_unlocked(mut __fp: *mut FILE) -> c_int {
    return if 0 !=
                  ((*__fp)._IO_read_ptr >= (*__fp)._IO_read_end) as
                      c_int as c_long {
               __uflow(__fp)
           } else {
               let fresh0 = (*__fp)._IO_read_ptr;
               (*__fp)._IO_read_ptr = (*__fp)._IO_read_ptr.offset(1);
               *(fresh0 as *mut c_uchar) as c_int
           };
}
unsafe extern "C" fn getchar_unlocked() -> c_int {
    return if 0 !=
                  ((*stdin)._IO_read_ptr >= (*stdin)._IO_read_end) as
                      c_int as c_long {
               __uflow(stdin)
           } else {
               let fresh1 = (*stdin)._IO_read_ptr;
               (*stdin)._IO_read_ptr = (*stdin)._IO_read_ptr.offset(1);
               *(fresh1 as *mut c_uchar) as c_int
           };
}
unsafe extern "C" fn fgetc_unlocked(mut __fp: *mut FILE) -> c_int {
    return if 0 !=
                  ((*__fp)._IO_read_ptr >= (*__fp)._IO_read_end) as
                      c_int as c_long {
               __uflow(__fp)
           } else {
               let fresh2 = (*__fp)._IO_read_ptr;
               (*__fp)._IO_read_ptr = (*__fp)._IO_read_ptr.offset(1);
               *(fresh2 as *mut c_uchar) as c_int
           };
}
unsafe extern "C" fn putchar(mut __c: c_int) -> c_int {
    return _IO_putc(__c, stdout);
}
unsafe extern "C" fn fputc_unlocked(mut __c: c_int,
                                    mut __stream: *mut FILE) -> c_int {
    return if 0 !=
                  ((*__stream)._IO_write_ptr >= (*__stream)._IO_write_end) as
                      c_int as c_long {
               __overflow(__stream, __c as c_uchar as c_int)
           } else {
               let fresh3 = (*__stream)._IO_write_ptr;
               (*__stream)._IO_write_ptr =
                   (*__stream)._IO_write_ptr.offset(1);
               *fresh3 = __c as c_char;
               *fresh3 as c_uchar as c_int
           };
}
unsafe extern "C" fn putc_unlocked(mut __c: c_int,
                                   mut __stream: *mut FILE) -> c_int {
    return if 0 !=
                  ((*__stream)._IO_write_ptr >= (*__stream)._IO_write_end) as
                      c_int as c_long {
               __overflow(__stream, __c as c_uchar as c_int)
           } else {
               let fresh4 = (*__stream)._IO_write_ptr;
               (*__stream)._IO_write_ptr =
                   (*__stream)._IO_write_ptr.offset(1);
               *fresh4 = __c as c_char;
               *fresh4 as c_uchar as c_int
           };
}
unsafe extern "C" fn putchar_unlocked(mut __c: c_int) -> c_int {
    return if 0 !=
                  ((*stdout)._IO_write_ptr >= (*stdout)._IO_write_end) as
                      c_int as c_long {
               __overflow(stdout, __c as c_uchar as c_int)
           } else {
               let fresh5 = (*stdout)._IO_write_ptr;
               (*stdout)._IO_write_ptr = (*stdout)._IO_write_ptr.offset(1);
               *fresh5 = __c as c_char;
               *fresh5 as c_uchar as c_int
           };
}
unsafe extern "C" fn feof_unlocked(mut __stream: *mut FILE) -> c_int {
    return ((*__stream)._flags & 16i32 != 0i32) as c_int;
}
unsafe extern "C" fn ferror_unlocked(mut __stream: *mut FILE) -> c_int {
    return ((*__stream)._flags & 32i32 != 0i32) as c_int;
}
unsafe extern "C" fn Rans64MulHi(mut a: uint64_t, mut b: uint64_t)
 -> uint64_t {
    return ((a as u128).wrapping_mul(b as u128) >> 64i32) as uint64_t;
}
unsafe extern "C" fn Rans64EncInit(mut r: *mut Rans64State) -> () {
    *r = (1u64 << 31i32) as Rans64State;
}
unsafe extern "C" fn Rans64EncPut(mut r: *mut Rans64State,
                                  mut pptr: *mut *mut uint32_t,
                                  mut start: uint32_t, mut freq: uint32_t,
                                  mut scale_bits: uint32_t) -> () {
    let mut x: uint64_t = *r;
    let mut x_max: uint64_t =
        (1u64 << 31i32 >> scale_bits <<
             32i32).wrapping_mul(freq as c_ulonglong) as uint64_t;
    if x >= x_max {
        *pptr = (*pptr).offset(-1isize);
        **pptr = x as uint32_t;
        x >>= 32i32
    }
    *r =
        (x.wrapping_div(freq as c_ulong) <<
             scale_bits).wrapping_add(x.wrapping_rem(freq as
                                                         c_ulong)).wrapping_add(start
                                                                                          as
                                                                                          c_ulong);
}
unsafe extern "C" fn Rans64EncFlush(mut r: *mut Rans64State,
                                    mut pptr: *mut *mut uint32_t) -> () {
    let mut x: uint64_t = *r;
    *pptr = (*pptr).offset(-2isize);
    *(*pptr).offset(0isize) = (x >> 0i32) as uint32_t;
    *(*pptr).offset(1isize) = (x >> 32i32) as uint32_t;
}
unsafe extern "C" fn Rans64DecInit(mut r: *mut Rans64State,
                                   mut pptr: *mut *mut uint32_t) -> () {
    let mut x: uint64_t = 0;
    x = (*(*pptr).offset(0isize) as uint64_t) << 0i32;
    x |= (*(*pptr).offset(1isize) as uint64_t) << 32i32;
    *pptr = (*pptr).offset(2isize);
    *r = x;
}
unsafe extern "C" fn Rans64DecGet(mut r: *mut Rans64State,
                                  mut scale_bits: uint32_t) -> uint32_t {
    return (*r &
                (1u32 << scale_bits).wrapping_sub(1i32 as c_uint) as
                    c_ulong) as uint32_t;
}
unsafe extern "C" fn Rans64DecAdvance(mut r: *mut Rans64State,
                                      mut pptr: *mut *mut uint32_t,
                                      mut start: uint32_t, mut freq: uint32_t,
                                      mut scale_bits: uint32_t) -> () {
    let mut mask: uint64_t =
        (1u64 << scale_bits).wrapping_sub(1i32 as c_ulonglong) as
            uint64_t;
    let mut x: uint64_t = *r;
    x =
        (freq as
             c_ulong).wrapping_mul(x >>
                                             scale_bits).wrapping_add(x &
                                                                          mask).wrapping_sub(start
                                                                                                 as
                                                                                                 c_ulong);
    if (x as c_ulonglong) < 1u64 << 31i32 {
        x = x << 32i32 | **pptr as c_ulong;
        *pptr = (*pptr).offset(1isize)
    }
    *r = x;
}
unsafe extern "C" fn Rans64EncSymbolInit(mut s: *mut Rans64EncSymbol,
                                         mut start: uint32_t,
                                         mut freq: uint32_t,
                                         mut scale_bits: uint32_t) -> () {
    (*s).freq = freq;
    (*s).cmpl_freq =
        ((1i32 << scale_bits) as c_uint).wrapping_sub(freq);
    if freq < 2i32 as c_uint {
        (*s).rcp_freq = !0u64 as uint64_t;
        (*s).rcp_shift = 0i32 as uint32_t;
        (*s).bias =
            start.wrapping_add((1i32 << scale_bits) as
                                   c_uint).wrapping_sub(1i32 as
                                                                  c_uint)
    } else {
        let mut shift: uint32_t = 0i32 as uint32_t;
        let mut x0: uint64_t = 0;
        let mut x1: uint64_t = 0;
        let mut t0: uint64_t = 0;
        let mut t1: uint64_t = 0;
        while freq > 1u32 << shift { shift = shift.wrapping_add(1) }
        x0 = freq.wrapping_sub(1i32 as c_uint) as uint64_t;
        x1 = (1u64 << shift.wrapping_add(31i32 as c_uint)) as uint64_t;
        t1 = x1.wrapping_div(freq as c_ulong);
        x0 =
            (x0 as
                 c_ulong).wrapping_add(x1.wrapping_rem(freq as
                                                                 c_ulong)
                                                 << 32i32) as uint64_t as
                uint64_t;
        t0 = x0.wrapping_div(freq as c_ulong);
        (*s).rcp_freq = t0.wrapping_add(t1 << 32i32);
        (*s).rcp_shift = shift.wrapping_sub(1i32 as c_uint);
        (*s).bias = start
    };
}
unsafe extern "C" fn Rans64DecSymbolInit(mut s: *mut Rans64DecSymbol,
                                         mut start: uint32_t,
                                         mut freq: uint32_t) -> () {
    (*s).start = start;
    (*s).freq = freq;
}
unsafe extern "C" fn Rans64EncPutSymbol(mut r: *mut Rans64State,
                                        mut pptr: *mut *mut uint32_t,
                                        mut sym: *const Rans64EncSymbol,
                                        mut scale_bits: uint32_t) -> () {
    let mut x: uint64_t = *r;
    let mut x_max: uint64_t =
        (1u64 << 31i32 >> scale_bits <<
             32i32).wrapping_mul((*sym).freq as c_ulonglong) as
            uint64_t;
    if x >= x_max {
        *pptr = (*pptr).offset(-1isize);
        **pptr = x as uint32_t;
        x >>= 32i32
    }
    let mut q: uint64_t = Rans64MulHi(x, (*sym).rcp_freq) >> (*sym).rcp_shift;
    *r =
        x.wrapping_add((*sym).bias as
                           c_ulong).wrapping_add(q.wrapping_mul((*sym).cmpl_freq
                                                                          as
                                                                          c_ulong));
}
unsafe extern "C" fn Rans64DecAdvanceSymbol(mut r: *mut Rans64State,
                                            mut pptr: *mut *mut uint32_t,
                                            mut sym: *const Rans64DecSymbol,
                                            mut scale_bits: uint32_t) -> () {
    Rans64DecAdvance(r, pptr, (*sym).start, (*sym).freq, scale_bits);
}
unsafe extern "C" fn Rans64DecAdvanceStep(mut r: *mut Rans64State,
                                          mut start: uint32_t,
                                          mut freq: uint32_t,
                                          mut scale_bits: uint32_t) -> () {
    let mut mask: uint64_t =
        (1u32 << scale_bits).wrapping_sub(1i32 as c_uint) as uint64_t;
    let mut x: uint64_t = *r;
    *r =
        (freq as
             c_ulong).wrapping_mul(x >>
                                             scale_bits).wrapping_add(x &
                                                                          mask).wrapping_sub(start
                                                                                                 as
                                                                                                 c_ulong);
}
unsafe extern "C" fn Rans64DecAdvanceSymbolStep(mut r: *mut Rans64State,
                                                mut sym:
                                                    *const Rans64DecSymbol,
                                                mut scale_bits: uint32_t)
 -> () {
    Rans64DecAdvanceStep(r, (*sym).start, (*sym).freq, scale_bits);
}
unsafe extern "C" fn Rans64DecRenorm(mut r: *mut Rans64State,
                                     mut pptr: *mut *mut uint32_t) -> () {
    let mut x: uint64_t = *r;
    if (x as c_ulonglong) < 1u64 << 31i32 {
        x = x << 32i32 | **pptr as c_ulong;
        *pptr = (*pptr).offset(1isize)
    }
    *r = x;
}
unsafe extern "C" fn atof(mut __nptr: *const c_char) -> c_double {
    return strtod(__nptr, 0 as *mut c_void as *mut *mut c_char);
}
unsafe extern "C" fn atoi(mut __nptr: *const c_char) -> c_int {
    return strtol(__nptr, 0 as *mut c_void as *mut *mut c_char,
                  10i32) as c_int;
}
unsafe extern "C" fn atol(mut __nptr: *const c_char) -> c_long {
    return strtol(__nptr, 0 as *mut c_void as *mut *mut c_char,
                  10i32);
}
unsafe extern "C" fn atoll(mut __nptr: *const c_char)
 -> c_longlong {
    return strtoll(__nptr, 0 as *mut c_void as *mut *mut c_char,
                   10i32);
}
unsafe extern "C" fn gnu_dev_major(mut __dev: c_ulonglong)
 -> c_uint {
    return (__dev >> 8i32 & 4095i32 as c_ulonglong |
                ((__dev >> 32i32) as c_uint & !4095i32 as c_uint)
                    as c_ulonglong) as c_uint;
}
unsafe extern "C" fn gnu_dev_minor(mut __dev: c_ulonglong)
 -> c_uint {
    return (__dev & 255i32 as c_ulonglong |
                ((__dev >> 12i32) as c_uint & !255i32 as c_uint)
                    as c_ulonglong) as c_uint;
}
unsafe extern "C" fn gnu_dev_makedev(mut __major: c_uint,
                                     mut __minor: c_uint)
 -> c_ulonglong {
    return (__minor & 255i32 as c_uint |
                (__major & 4095i32 as c_uint) << 8i32) as
               c_ulonglong |
               ((__minor & !255i32 as c_uint) as c_ulonglong) <<
                   12i32 |
               ((__major & !4095i32 as c_uint) as c_ulonglong) <<
                   32i32;
}
unsafe extern "C" fn bsearch(mut __key: *const c_void,
                             mut __base: *const c_void,
                             mut __nmemb: size_t, mut __size: size_t,
                             mut __compar: __compar_fn_t)
 -> *mut c_void {
    let mut __l: size_t = 0;
    let mut __u: size_t = 0;
    let mut __idx: size_t = 0;
    let mut __p: *const c_void = 0 as *const c_void;
    let mut __comparison: c_int = 0;
    __l = 0i32 as size_t;
    __u = __nmemb;
    while __l < __u {
        __idx = __l.wrapping_add(__u).wrapping_div(2i32 as c_ulong);
        __p =
            (__base as
                 *const c_char).offset(__idx.wrapping_mul(__size) as
                                                 isize) as *mut c_void;
        __comparison =
            __compar.expect("non-null function pointer")(__key, __p);
        if __comparison < 0i32 {
            __u = __idx
        } else if __comparison > 0i32 {
            __l = __idx.wrapping_add(1i32 as c_ulong)
        } else { return __p as *mut c_void }
    }
    return 0 as *mut c_void;
}
unsafe extern "C" fn __strcspn_c1(mut __s: *const c_char,
                                  mut __reject: c_int) -> size_t {
    let mut __result: size_t = 0i32 as size_t;
    while *__s.offset(__result as isize) as c_int != '\u{0}' as i32 &&
              *__s.offset(__result as isize) as c_int != __reject {
        __result = __result.wrapping_add(1)
    }
    return __result;
}
unsafe extern "C" fn __strcspn_c2(mut __s: *const c_char,
                                  mut __reject1: c_int,
                                  mut __reject2: c_int) -> size_t {
    let mut __result: size_t = 0i32 as size_t;
    while *__s.offset(__result as isize) as c_int != '\u{0}' as i32 &&
              *__s.offset(__result as isize) as c_int != __reject1 &&
              *__s.offset(__result as isize) as c_int != __reject2 {
        __result = __result.wrapping_add(1)
    }
    return __result;
}
unsafe extern "C" fn __strcspn_c3(mut __s: *const c_char,
                                  mut __reject1: c_int,
                                  mut __reject2: c_int,
                                  mut __reject3: c_int) -> size_t {
    let mut __result: size_t = 0i32 as size_t;
    while *__s.offset(__result as isize) as c_int != '\u{0}' as i32 &&
              *__s.offset(__result as isize) as c_int != __reject1 &&
              *__s.offset(__result as isize) as c_int != __reject2 &&
              *__s.offset(__result as isize) as c_int != __reject3 {
        __result = __result.wrapping_add(1)
    }
    return __result;
}
unsafe extern "C" fn __strspn_c1(mut __s: *const c_char,
                                 mut __accept: c_int) -> size_t {
    let mut __result: size_t = 0i32 as size_t;
    while *__s.offset(__result as isize) as c_int == __accept {
        __result = __result.wrapping_add(1)
    }
    return __result;
}
unsafe extern "C" fn __strspn_c2(mut __s: *const c_char,
                                 mut __accept1: c_int,
                                 mut __accept2: c_int) -> size_t {
    let mut __result: size_t = 0i32 as size_t;
    while *__s.offset(__result as isize) as c_int == __accept1 ||
              *__s.offset(__result as isize) as c_int == __accept2 {
        __result = __result.wrapping_add(1)
    }
    return __result;
}
unsafe extern "C" fn __strspn_c3(mut __s: *const c_char,
                                 mut __accept1: c_int,
                                 mut __accept2: c_int,
                                 mut __accept3: c_int) -> size_t {
    let mut __result: size_t = 0i32 as size_t;
    while *__s.offset(__result as isize) as c_int == __accept1 ||
              *__s.offset(__result as isize) as c_int == __accept2 ||
              *__s.offset(__result as isize) as c_int == __accept3 {
        __result = __result.wrapping_add(1)
    }
    return __result;
}
unsafe extern "C" fn __strpbrk_c2(mut __s: *const c_char,
                                  mut __accept1: c_int,
                                  mut __accept2: c_int)
 -> *mut c_char {
    while *__s as c_int != '\u{0}' as i32 &&
              *__s as c_int != __accept1 &&
              *__s as c_int != __accept2 {
        __s = __s.offset(1isize)
    }
    return if *__s as c_int == '\u{0}' as i32 {
               0 as *mut c_char
           } else { __s as size_t as *mut c_char };
}
unsafe extern "C" fn __strpbrk_c3(mut __s: *const c_char,
                                  mut __accept1: c_int,
                                  mut __accept2: c_int,
                                  mut __accept3: c_int)
 -> *mut c_char {
    while *__s as c_int != '\u{0}' as i32 &&
              *__s as c_int != __accept1 &&
              *__s as c_int != __accept2 &&
              *__s as c_int != __accept3 {
        __s = __s.offset(1isize)
    }
    return if *__s as c_int == '\u{0}' as i32 {
               0 as *mut c_char
           } else { __s as size_t as *mut c_char };
}
unsafe extern "C" fn __strtok_r_1c(mut __s: *mut c_char,
                                   mut __sep: c_char,
                                   mut __nextp: *mut *mut c_char)
 -> *mut c_char {
    let mut __result: *mut c_char = 0 as *mut c_char;
    if __s.is_null() { __s = *__nextp }
    while *__s as c_int == __sep as c_int {
        __s = __s.offset(1isize)
    }
    __result = 0 as *mut c_char;
    if *__s as c_int != '\u{0}' as i32 {
        let fresh6 = __s;
        __s = __s.offset(1);
        __result = fresh6;
        while *__s as c_int != '\u{0}' as i32 {
            let fresh7 = __s;
            __s = __s.offset(1);
            if !(*fresh7 as c_int == __sep as c_int) {
                continue ;
            }
            *__s.offset(-1i32 as isize) = '\u{0}' as i32 as c_char;
            break ;
        }
    }
    *__nextp = __s;
    return __result;
}
unsafe extern "C" fn __strsep_2c(mut __s: *mut *mut c_char,
                                 mut __reject1: c_char,
                                 mut __reject2: c_char)
 -> *mut c_char {
    let mut __retval: *mut c_char = *__s;
    if !__retval.is_null() {
        let mut __cp: *mut c_char = __retval;
        loop  {
            if *__cp as c_int == '\u{0}' as i32 {
                __cp = 0 as *mut c_char;
                break ;
            } else if *__cp as c_int == __reject1 as c_int ||
                          *__cp as c_int == __reject2 as c_int {
                let fresh8 = __cp;
                __cp = __cp.offset(1);
                *fresh8 = '\u{0}' as i32 as c_char;
                break ;
            } else { __cp = __cp.offset(1isize) }
        }
        *__s = __cp
    }
    return __retval;
}
unsafe extern "C" fn __strsep_3c(mut __s: *mut *mut c_char,
                                 mut __reject1: c_char,
                                 mut __reject2: c_char,
                                 mut __reject3: c_char)
 -> *mut c_char {
    let mut __retval: *mut c_char = *__s;
    if !__retval.is_null() {
        let mut __cp: *mut c_char = __retval;
        loop  {
            if *__cp as c_int == '\u{0}' as i32 {
                __cp = 0 as *mut c_char;
                break ;
            } else if *__cp as c_int == __reject1 as c_int ||
                          *__cp as c_int == __reject2 as c_int ||
                          *__cp as c_int == __reject3 as c_int {
                let fresh9 = __cp;
                __cp = __cp.offset(1);
                *fresh9 = '\u{0}' as i32 as c_char;
                break ;
            } else { __cp = __cp.offset(1isize) }
        }
        *__s = __cp
    }
    return __retval;
}
unsafe extern "C" fn hist8(mut in_0: *mut c_uchar,
                           mut in_size: c_uint,
                           mut F0: *mut c_int) -> () {
    let mut F1: [c_int; 264] =
        [0i32, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
    let mut F2: [c_int; 264] =
        [0i32, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
    let mut F3: [c_int; 264] =
        [0i32, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
    let mut F4: [c_int; 264] =
        [0i32, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
    let mut F5: [c_int; 264] =
        [0i32, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
    let mut F6: [c_int; 264] =
        [0i32, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
    let mut F7: [c_int; 264] =
        [0i32, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
    let mut i: c_int = 0;
    let mut i8: c_int =
        (in_size & !8i32 as c_uint) as c_int;
    i = 0i32;
    while i < i8 {
        let ref mut fresh10 =
            *F0.offset(*in_0.offset((i + 0i32) as isize) as isize);
        *fresh10 += 1;
        F1[*in_0.offset((i + 1i32) as isize) as usize] += 1;
        F2[*in_0.offset((i + 2i32) as isize) as usize] += 1;
        F3[*in_0.offset((i + 3i32) as isize) as usize] += 1;
        F4[*in_0.offset((i + 4i32) as isize) as usize] += 1;
        F5[*in_0.offset((i + 5i32) as isize) as usize] += 1;
        F6[*in_0.offset((i + 6i32) as isize) as usize] += 1;
        F7[*in_0.offset((i + 7i32) as isize) as usize] += 1;
        i += 8i32
    }
    while (i as c_uint) < in_size {
        let fresh11 = i;
        i = i + 1;
        let ref mut fresh12 =
            *F0.offset(*in_0.offset(fresh11 as isize) as isize);
        *fresh12 += 1
    }
    i = 0i32;
    while i < 256i32 {
        *F0.offset(i as isize) +=
            F1[i as usize] + F2[i as usize] + F3[i as usize] + F4[i as usize]
                + F5[i as usize] + F6[i as usize] + F7[i as usize];
        i += 1
    };
}
#[no_mangle]
pub unsafe extern "C" fn rans_compress_O0(mut in_0: *mut c_uchar,
                                          mut in_size: c_uint,
                                          mut out_size: *mut c_uint)
 -> *mut c_uchar {
    let mut current_block: u64;
    let mut sz: c_int =
        (1.05f64 * in_size as c_double +
             (257i32 * 257i32 * 3i32) as c_double +
             8i32 as c_double) as c_int & !3i32;
    let mut out_buf: *mut c_uchar =
        malloc((sz + 4i32) as c_ulong) as *mut c_uchar;
    let mut cp: *mut c_uchar = 0 as *mut c_uchar;
    let mut out_end: *mut c_uchar = 0 as *mut c_uchar;
    let mut syms: [RansEncSymbol; 256] =
        [RansEncSymbol{rcp_freq: 0,
                       freq: 0,
                       bias: 0,
                       cmpl_freq: 0,
                       rcp_shift: 0,}; 256];
    let mut rans0: RansState = 0;
    let mut rans1: RansState = 0;
    let mut rans2: RansState = 0;
    let mut rans3: RansState = 0;
    let mut ptr: *mut uint8_t = 0 as *mut uint8_t;
    let mut F: [c_int; 264] =
        [0i32, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
    let mut i: c_int = 0;
    let mut j: c_int = 0;
    let mut tab_size: c_int = 0;
    let mut rle: c_int = 0;
    let mut x: c_int = 0;
    let mut fsum: c_int = 0i32;
    let mut m: c_int = 0i32;
    let mut M: c_int = 0i32;
    let mut tr: uint64_t = 0;
    if out_buf.is_null() {
        return 0 as *mut c_uchar
    } else {
        out_end = out_buf.offset(sz as isize);
        ptr = out_end;
        hist8(in_0, in_size, F.as_mut_ptr());
        tr =
            (((1i32 << 12i32) as uint64_t) <<
                 31i32).wrapping_div(in_size as
                                         c_ulong).wrapping_add(((1i32 <<
                                                                           30i32)
                                                                          as
                                                                          c_uint).wrapping_div(in_size)
                                                                         as
                                                                         c_ulong);
        j = 0i32;
        M = j;
        m = M;
        while j < 256i32 {
            if !(0 == F[j as usize]) {
                if m < F[j as usize] { m = F[j as usize]; M = j }
                F[j as usize] =
                    ((F[j as usize] as c_ulong).wrapping_mul(tr) >>
                         31i32) as c_int;
                if F[j as usize] == 0i32 { F[j as usize] = 1i32 }
                fsum += F[j as usize]
            }
            j += 1
        }
        fsum += 1;
        if fsum < 1i32 << 12i32 {
            F[M as usize] += (1i32 << 12i32) - fsum
        } else { F[M as usize] -= fsum - (1i32 << 12i32) }
        cp = out_buf.offset(4isize);
        j = 0i32;
        rle = j;
        x = rle;
        while j < 256i32 {
            if 0 != F[j as usize] {
                if 0 != rle {
                    rle -= 1
                } else {
                    let fresh13 = cp;
                    cp = cp.offset(1);
                    *fresh13 = j as c_uchar;
                    if 0 == rle && 0 != j && 0 != F[(j - 1i32) as usize] {
                        rle = j + 1i32;
                        while rle < 256i32 && 0 != F[rle as usize] {
                            rle += 1
                        }
                        rle -= j + 1i32;
                        let fresh14 = cp;
                        cp = cp.offset(1);
                        *fresh14 = rle as c_uchar
                    }
                }
                if F[j as usize] < 128i32 {
                    let fresh15 = cp;
                    cp = cp.offset(1);
                    *fresh15 = F[j as usize] as c_uchar
                } else {
                    let fresh16 = cp;
                    cp = cp.offset(1);
                    *fresh16 =
                        (128i32 | F[j as usize] >> 8i32) as c_uchar;
                    let fresh17 = cp;
                    cp = cp.offset(1);
                    *fresh17 = (F[j as usize] & 255i32) as c_uchar
                }
                Rans64EncSymbolInit(&mut syms[j as usize] as
                                        *mut RansEncSymbol, x as uint32_t,
                                    F[j as usize] as uint32_t,
                                    12i32 as uint32_t);
                x += F[j as usize]
            }
            j += 1
        }
        let fresh18 = cp;
        cp = cp.offset(1);
        *fresh18 = 0i32 as c_uchar;
        tab_size =
            out_buf.offset_to(cp).expect("bad offset_to") as c_long as
                c_int;
        Rans64EncInit(&mut rans0 as *mut RansState);
        Rans64EncInit(&mut rans1 as *mut RansState);
        Rans64EncInit(&mut rans2 as *mut RansState);
        Rans64EncInit(&mut rans3 as *mut RansState);
        i = (in_size & 3i32 as c_uint) as c_int;
        match i {
            3 => {
                Rans64EncPutSymbol(&mut rans2 as *mut RansState,
                                   &mut ptr as *mut *mut uint8_t as
                                       *mut *mut uint32_t,
                                   &mut syms[*in_0.offset(in_size.wrapping_sub((i
                                                                                    -
                                                                                    2i32)
                                                                                   as
                                                                                   c_uint)
                                                              as isize) as
                                                 usize] as *mut RansEncSymbol,
                                   12i32 as uint32_t);
                current_block = 14453171571987214953;
            }
            2 => { current_block = 14453171571987214953; }
            1 => { current_block = 15409834050583242653; }
            0 | _ => { current_block = 224731115979188411; }
        }
        match current_block {
            14453171571987214953 => {
                Rans64EncPutSymbol(&mut rans1 as *mut RansState,
                                   &mut ptr as *mut *mut uint8_t as
                                       *mut *mut uint32_t,
                                   &mut syms[*in_0.offset(in_size.wrapping_sub((i
                                                                                    -
                                                                                    1i32)
                                                                                   as
                                                                                   c_uint)
                                                              as isize) as
                                                 usize] as *mut RansEncSymbol,
                                   12i32 as uint32_t);
                current_block = 15409834050583242653;
            }
            _ => { }
        }
        match current_block {
            15409834050583242653 => {
                Rans64EncPutSymbol(&mut rans0 as *mut RansState,
                                   &mut ptr as *mut *mut uint8_t as
                                       *mut *mut uint32_t,
                                   &mut syms[*in_0.offset(in_size.wrapping_sub((i
                                                                                    -
                                                                                    0i32)
                                                                                   as
                                                                                   c_uint)
                                                              as isize) as
                                                 usize] as *mut RansEncSymbol,
                                   12i32 as uint32_t);
            }
            _ => { }
        }
        i = (in_size & !3i32 as c_uint) as c_int;
        while i > 0i32 {
            let mut s3: *mut RansEncSymbol =
                &mut syms[*in_0.offset((i - 1i32) as isize) as usize] as
                    *mut RansEncSymbol;
            let mut s2: *mut RansEncSymbol =
                &mut syms[*in_0.offset((i - 2i32) as isize) as usize] as
                    *mut RansEncSymbol;
            let mut s1: *mut RansEncSymbol =
                &mut syms[*in_0.offset((i - 3i32) as isize) as usize] as
                    *mut RansEncSymbol;
            let mut s0: *mut RansEncSymbol =
                &mut syms[*in_0.offset((i - 4i32) as isize) as usize] as
                    *mut RansEncSymbol;
            Rans64EncPutSymbol(&mut rans3 as *mut RansState,
                               &mut ptr as *mut *mut uint8_t as
                                   *mut *mut uint32_t, s3, 12i32 as uint32_t);
            Rans64EncPutSymbol(&mut rans2 as *mut RansState,
                               &mut ptr as *mut *mut uint8_t as
                                   *mut *mut uint32_t, s2, 12i32 as uint32_t);
            Rans64EncPutSymbol(&mut rans1 as *mut RansState,
                               &mut ptr as *mut *mut uint8_t as
                                   *mut *mut uint32_t, s1, 12i32 as uint32_t);
            Rans64EncPutSymbol(&mut rans0 as *mut RansState,
                               &mut ptr as *mut *mut uint8_t as
                                   *mut *mut uint32_t, s0, 12i32 as uint32_t);
            i -= 4i32
        }
        Rans64EncFlush(&mut rans3 as *mut RansState,
                       &mut ptr as *mut *mut uint8_t as *mut *mut uint32_t);
        Rans64EncFlush(&mut rans2 as *mut RansState,
                       &mut ptr as *mut *mut uint8_t as *mut *mut uint32_t);
        Rans64EncFlush(&mut rans1 as *mut RansState,
                       &mut ptr as *mut *mut uint8_t as *mut *mut uint32_t);
        Rans64EncFlush(&mut rans0 as *mut RansState,
                       &mut ptr as *mut *mut uint8_t as *mut *mut uint32_t);
        *out_size =
            (ptr.offset_to(out_end).expect("bad offset_to") as c_long +
                 tab_size as c_long) as c_uint;
        cp = out_buf;
        let fresh19 = cp;
        cp = cp.offset(1);
        *fresh19 =
            (in_size >> 0i32 & 255i32 as c_uint) as c_uchar;
        let fresh20 = cp;
        cp = cp.offset(1);
        *fresh20 =
            (in_size >> 8i32 & 255i32 as c_uint) as c_uchar;
        let fresh21 = cp;
        cp = cp.offset(1);
        *fresh21 =
            (in_size >> 16i32 & 255i32 as c_uint) as c_uchar;
        let fresh22 = cp;
        cp = cp.offset(1);
        *fresh22 =
            (in_size >> 24i32 & 255i32 as c_uint) as c_uchar;
        memmove(out_buf.offset(tab_size as isize) as *mut c_void,
                ptr as *const c_void,
                ptr.offset_to(out_end).expect("bad offset_to") as c_long
                    as c_ulong);
        return out_buf
    };
}
#[no_mangle]
pub unsafe extern "C" fn rans_uncompress_O0(mut in_0: *mut c_uchar,
                                            mut in_size: c_uint,
                                            mut out_size: *mut c_uint)
 -> *mut c_uchar {
    let mut c_0: c_uchar = 0;
    let mut cp: *mut c_uchar = in_0.offset(4isize);
    let mut i: c_int = 0;
    let mut j: c_int = 0;
    let mut x: c_int = 0;
    let mut out_sz: c_int = 0;
    let mut rle: c_int = 0;
    let mut out_buf: *mut c_char = 0 as *mut c_char;
    let mut D: ari_decoder = ari_decoder{R: [0; 4096],};
    let mut syms: [RansDecSymbol; 256] =
        [RansDecSymbol{start: 0, freq: 0,}; 256];
    out_sz =
        (*in_0.offset(0isize) as c_int) << 0i32 |
            (*in_0.offset(1isize) as c_int) << 8i32 |
            (*in_0.offset(2isize) as c_int) << 16i32 |
            (*in_0.offset(3isize) as c_int) << 24i32;
    out_buf = malloc(out_sz as c_ulong) as *mut c_char;
    if out_buf.is_null() {
        return 0 as *mut c_uchar
    } else {
        x = 0i32;
        rle = x;
        let fresh23 = cp;
        cp = cp.offset(1);
        j = *fresh23 as c_int;
        loop  {
            let mut F: c_int = 0;
            let mut C: c_int = 0;
            let fresh24 = cp;
            cp = cp.offset(1);
            F = *fresh24 as c_int;
            if F >= 128i32 {
                F &= !128i32;
                let fresh25 = cp;
                cp = cp.offset(1);
                F = (F & 127i32) << 8i32 | *fresh25 as c_int
            }
            C = x;
            Rans64DecSymbolInit(&mut syms[j as usize] as *mut RansDecSymbol,
                                C as uint32_t, F as uint32_t);
            memset(&mut D.R[x as usize] as *mut c_uchar as
                       *mut c_void, j, F as c_ulong);
            x += F;
            if 0 == rle && j + 1i32 == *cp as c_int {
                let fresh26 = cp;
                cp = cp.offset(1);
                j = *fresh26 as c_int;
                let fresh27 = cp;
                cp = cp.offset(1);
                rle = *fresh27 as c_int
            } else if 0 != rle {
                rle -= 1;
                j += 1
            } else {
                let fresh28 = cp;
                cp = cp.offset(1);
                j = *fresh28 as c_int
            }
            if !(0 != j) { break ; }
        }
        let mut rans0: RansState = 0;
        let mut rans1: RansState = 0;
        let mut rans2: RansState = 0;
        let mut rans3: RansState = 0;
        let mut ptr: *mut uint8_t = cp;
        Rans64DecInit(&mut rans0 as *mut RansState,
                      &mut ptr as *mut *mut uint8_t as *mut *mut uint32_t);
        Rans64DecInit(&mut rans1 as *mut RansState,
                      &mut ptr as *mut *mut uint8_t as *mut *mut uint32_t);
        Rans64DecInit(&mut rans2 as *mut RansState,
                      &mut ptr as *mut *mut uint8_t as *mut *mut uint32_t);
        Rans64DecInit(&mut rans3 as *mut RansState,
                      &mut ptr as *mut *mut uint8_t as *mut *mut uint32_t);
        let mut out_end: c_int = out_sz & !3i32;
        let mut R: [RansState; 4] = [0; 4];
        R[0usize] = rans0;
        R[1usize] = rans1;
        R[2usize] = rans2;
        R[3usize] = rans3;
        let mut mask: uint32_t =
            (1u32 << 12i32).wrapping_sub(1i32 as c_uint);
        i = 0i32;
        while i < out_end {
            let mut m: [uint32_t; 4] = [0; 4];
            let mut c: [uint8_t; 4] = [0; 4];
            m[0usize] = (R[0usize] & mask as c_ulong) as uint32_t;
            c[0usize] = D.R[m[0usize] as usize];
            *out_buf.offset((i + 0i32) as isize) = c[0usize] as c_char;
            R[0usize] =
                (syms[c[0usize] as usize].freq as
                     c_ulong).wrapping_mul(R[0usize] >> 12i32);
            R[0usize] =
                (R[0usize] as
                     c_ulong).wrapping_add(m[0usize].wrapping_sub(syms[c[0usize]
                                                                                 as
                                                                                 usize].start)
                                                     as c_ulong) as
                    RansState as RansState;
            m[1usize] = (R[1usize] & mask as c_ulong) as uint32_t;
            c[1usize] = D.R[m[1usize] as usize];
            *out_buf.offset((i + 1i32) as isize) = c[1usize] as c_char;
            R[1usize] =
                (syms[c[1usize] as usize].freq as
                     c_ulong).wrapping_mul(R[1usize] >> 12i32);
            R[1usize] =
                (R[1usize] as
                     c_ulong).wrapping_add(m[1usize].wrapping_sub(syms[c[1usize]
                                                                                 as
                                                                                 usize].start)
                                                     as c_ulong) as
                    RansState as RansState;
            m[2usize] = (R[2usize] & mask as c_ulong) as uint32_t;
            c[2usize] = D.R[m[2usize] as usize];
            *out_buf.offset((i + 2i32) as isize) = c[2usize] as c_char;
            R[2usize] =
                (syms[c[2usize] as usize].freq as
                     c_ulong).wrapping_mul(R[2usize] >> 12i32);
            R[2usize] =
                (R[2usize] as
                     c_ulong).wrapping_add(m[2usize].wrapping_sub(syms[c[2usize]
                                                                                 as
                                                                                 usize].start)
                                                     as c_ulong) as
                    RansState as RansState;
            m[3usize] = (R[3usize] & mask as c_ulong) as uint32_t;
            c[3usize] = D.R[m[3usize] as usize];
            *out_buf.offset((i + 3i32) as isize) = c[3usize] as c_char;
            R[3usize] =
                (syms[c[3usize] as usize].freq as
                     c_ulong).wrapping_mul(R[3usize] >> 12i32);
            R[3usize] =
                (R[3usize] as
                     c_ulong).wrapping_add(m[3usize].wrapping_sub(syms[c[3usize]
                                                                                 as
                                                                                 usize].start)
                                                     as c_ulong) as
                    RansState as RansState;
            Rans64DecRenorm(&mut R[0usize] as *mut RansState,
                            &mut ptr as *mut *mut uint8_t as
                                *mut *mut uint32_t);
            Rans64DecRenorm(&mut R[1usize] as *mut RansState,
                            &mut ptr as *mut *mut uint8_t as
                                *mut *mut uint32_t);
            Rans64DecRenorm(&mut R[2usize] as *mut RansState,
                            &mut ptr as *mut *mut uint8_t as
                                *mut *mut uint32_t);
            Rans64DecRenorm(&mut R[3usize] as *mut RansState,
                            &mut ptr as *mut *mut uint8_t as
                                *mut *mut uint32_t);
            i += 4i32
        }
        rans0 = R[0usize];
        rans1 = R[1usize];
        rans2 = R[2usize];
        rans3 = R[3usize];
        match out_sz & 3i32 {
            1 => {
                c_0 =
                    D.R[Rans64DecGet(&mut rans0 as *mut RansState,
                                     12i32 as uint32_t) as usize];
                Rans64DecAdvanceSymbol(&mut rans0 as *mut RansState,
                                       &mut ptr as *mut *mut uint8_t as
                                           *mut *mut uint32_t,
                                       &mut syms[c_0 as usize] as
                                           *mut RansDecSymbol,
                                       12i32 as uint32_t);
                *out_buf.offset(out_end as isize) = c_0 as c_char
            }
            2 => {
                c_0 =
                    D.R[Rans64DecGet(&mut rans0 as *mut RansState,
                                     12i32 as uint32_t) as usize];
                Rans64DecAdvanceSymbol(&mut rans0 as *mut RansState,
                                       &mut ptr as *mut *mut uint8_t as
                                           *mut *mut uint32_t,
                                       &mut syms[c_0 as usize] as
                                           *mut RansDecSymbol,
                                       12i32 as uint32_t);
                *out_buf.offset(out_end as isize) = c_0 as c_char;
                c_0 =
                    D.R[Rans64DecGet(&mut rans1 as *mut RansState,
                                     12i32 as uint32_t) as usize];
                Rans64DecAdvanceSymbol(&mut rans1 as *mut RansState,
                                       &mut ptr as *mut *mut uint8_t as
                                           *mut *mut uint32_t,
                                       &mut syms[c_0 as usize] as
                                           *mut RansDecSymbol,
                                       12i32 as uint32_t);
                *out_buf.offset((out_end + 1i32) as isize) =
                    c_0 as c_char
            }
            3 => {
                c_0 =
                    D.R[Rans64DecGet(&mut rans0 as *mut RansState,
                                     12i32 as uint32_t) as usize];
                Rans64DecAdvanceSymbol(&mut rans0 as *mut RansState,
                                       &mut ptr as *mut *mut uint8_t as
                                           *mut *mut uint32_t,
                                       &mut syms[c_0 as usize] as
                                           *mut RansDecSymbol,
                                       12i32 as uint32_t);
                *out_buf.offset(out_end as isize) = c_0 as c_char;
                c_0 =
                    D.R[Rans64DecGet(&mut rans1 as *mut RansState,
                                     12i32 as uint32_t) as usize];
                Rans64DecAdvanceSymbol(&mut rans1 as *mut RansState,
                                       &mut ptr as *mut *mut uint8_t as
                                           *mut *mut uint32_t,
                                       &mut syms[c_0 as usize] as
                                           *mut RansDecSymbol,
                                       12i32 as uint32_t);
                *out_buf.offset((out_end + 1i32) as isize) =
                    c_0 as c_char;
                c_0 =
                    D.R[Rans64DecGet(&mut rans2 as *mut RansState,
                                     12i32 as uint32_t) as usize];
                Rans64DecAdvanceSymbol(&mut rans2 as *mut RansState,
                                       &mut ptr as *mut *mut uint8_t as
                                           *mut *mut uint32_t,
                                       &mut syms[c_0 as usize] as
                                           *mut RansDecSymbol,
                                       12i32 as uint32_t);
                *out_buf.offset((out_end + 2i32) as isize) =
                    c_0 as c_char
            }
            0 | _ => { }
        }
        *out_size = out_sz as c_uint;
        return out_buf as *mut c_uchar
    };
}
#[no_mangle]
pub unsafe extern "C" fn rans_uncompress_O0b(mut in_0: *mut c_uchar,
                                             mut in_size: c_uint,
                                             mut out_size: *mut c_uint)
 -> *mut c_uchar {
    let mut c_0: c_uchar = 0;
    let mut cp: *mut c_uchar = in_0.offset(4isize);
    let mut i: c_int = 0;
    let mut j: c_int = 0;
    let mut x: c_int = 0;
    let mut out_sz: c_int = 0;
    let mut rle: c_int = 0;
    let mut out_buf: *mut c_char = 0 as *mut c_char;
    let mut D: ari_decoder = ari_decoder{R: [0; 4096],};
    let mut syms: [RansDecSymbol; 256] =
        [RansDecSymbol{start: 0, freq: 0,}; 256];
    out_sz =
        (*in_0.offset(0isize) as c_int) << 0i32 |
            (*in_0.offset(1isize) as c_int) << 8i32 |
            (*in_0.offset(2isize) as c_int) << 16i32 |
            (*in_0.offset(3isize) as c_int) << 24i32;
    out_buf = malloc(out_sz as c_ulong) as *mut c_char;
    if out_buf.is_null() {
        return 0 as *mut c_uchar
    } else {
        x = 0i32;
        rle = x;
        let fresh29 = cp;
        cp = cp.offset(1);
        j = *fresh29 as c_int;
        loop  {
            let mut F: c_int = 0;
            let mut C: c_int = 0;
            let fresh30 = cp;
            cp = cp.offset(1);
            F = *fresh30 as c_int;
            if F >= 128i32 {
                F &= !128i32;
                let fresh31 = cp;
                cp = cp.offset(1);
                F = (F & 127i32) << 8i32 | *fresh31 as c_int
            }
            C = x;
            Rans64DecSymbolInit(&mut syms[j as usize] as *mut RansDecSymbol,
                                C as uint32_t, F as uint32_t);
            memset(&mut D.R[x as usize] as *mut c_uchar as
                       *mut c_void, j, F as c_ulong);
            x += F;
            if 0 == rle && j + 1i32 == *cp as c_int {
                let fresh33 = cp;
                cp = cp.offset(1);
                j = *fresh33 as c_int;
                let fresh34 = cp;
                cp = cp.offset(1);
                rle = *fresh34 as c_int
            } else if 0 != rle {
                rle -= 1;
                j += 1
            } else {
                let fresh35 = cp;
                cp = cp.offset(1);
                j = *fresh35 as c_int
            }
            if !(0 != j) { break ; }
        }
        let mut rans0: RansState = 0;
        let mut rans1: RansState = 0;
        let mut rans2: RansState = 0;
        let mut rans3: RansState = 0;
        let mut ptr: *mut uint8_t = cp;
        Rans64DecInit(&mut rans0 as *mut RansState,
                      &mut ptr as *mut *mut uint8_t as *mut *mut uint32_t);
        Rans64DecInit(&mut rans1 as *mut RansState,
                      &mut ptr as *mut *mut uint8_t as *mut *mut uint32_t);
        Rans64DecInit(&mut rans2 as *mut RansState,
                      &mut ptr as *mut *mut uint8_t as *mut *mut uint32_t);
        Rans64DecInit(&mut rans3 as *mut RansState,
                      &mut ptr as *mut *mut uint8_t as *mut *mut uint32_t);
        let mut out_end: c_int = out_sz & !3i32;
        let mut R: [RansState; 4] = [0; 4];
        R[0usize] = rans0;
        R[1usize] = rans1;
        R[2usize] = rans2;
        R[3usize] = rans3;
        let mut mask: uint32_t =
            (1u32 << 12i32).wrapping_sub(1i32 as c_uint);
        i = 0i32;
        while i < out_end {
            let mut m: [uint32_t; 4] =
                [(R[0usize] & mask as c_ulong) as uint32_t,
                 (R[1usize] & mask as c_ulong) as uint32_t,
                 (R[2usize] & mask as c_ulong) as uint32_t,
                 (R[3usize] & mask as c_ulong) as uint32_t];
            let mut c: [uint8_t; 4] =
                [D.R[m[0usize] as usize], D.R[m[1usize] as usize],
                 D.R[m[2usize] as usize], D.R[m[3usize] as usize]];
            *out_buf.offset((i + 0i32) as isize) = c[0usize] as c_char;
            *out_buf.offset((i + 1i32) as isize) = c[1usize] as c_char;
            *out_buf.offset((i + 2i32) as isize) = c[2usize] as c_char;
            *out_buf.offset((i + 3i32) as isize) = c[3usize] as c_char;
            R[0usize] =
                (syms[c[0usize] as usize].freq as
                     c_ulong).wrapping_mul(R[0usize] >> 12i32);
            R[0usize] =
                (R[0usize] as
                     c_ulong).wrapping_add(m[0usize].wrapping_sub(syms[c[0usize]
                                                                                 as
                                                                                 usize].start)
                                                     as c_ulong) as
                    RansState as RansState;
            R[1usize] =
                (syms[c[1usize] as usize].freq as
                     c_ulong).wrapping_mul(R[1usize] >> 12i32);
            R[1usize] =
                (R[1usize] as
                     c_ulong).wrapping_add(m[1usize].wrapping_sub(syms[c[1usize]
                                                                                 as
                                                                                 usize].start)
                                                     as c_ulong) as
                    RansState as RansState;
            R[2usize] =
                (syms[c[2usize] as usize].freq as
                     c_ulong).wrapping_mul(R[2usize] >> 12i32);
            R[2usize] =
                (R[2usize] as
                     c_ulong).wrapping_add(m[2usize].wrapping_sub(syms[c[2usize]
                                                                                 as
                                                                                 usize].start)
                                                     as c_ulong) as
                    RansState as RansState;
            R[3usize] =
                (syms[c[3usize] as usize].freq as
                     c_ulong).wrapping_mul(R[3usize] >> 12i32);
            R[3usize] =
                (R[3usize] as
                     c_ulong).wrapping_add(m[3usize].wrapping_sub(syms[c[3usize]
                                                                                 as
                                                                                 usize].start)
                                                     as c_ulong) as
                    RansState as RansState;
            Rans64DecRenorm(&mut R[0usize] as *mut RansState,
                            &mut ptr as *mut *mut uint8_t as
                                *mut *mut uint32_t);
            Rans64DecRenorm(&mut R[1usize] as *mut RansState,
                            &mut ptr as *mut *mut uint8_t as
                                *mut *mut uint32_t);
            Rans64DecRenorm(&mut R[2usize] as *mut RansState,
                            &mut ptr as *mut *mut uint8_t as
                                *mut *mut uint32_t);
            Rans64DecRenorm(&mut R[3usize] as *mut RansState,
                            &mut ptr as *mut *mut uint8_t as
                                *mut *mut uint32_t);
            i += 4i32
        }
        rans0 = R[0usize];
        rans1 = R[1usize];
        rans2 = R[2usize];
        rans3 = R[3usize];
        match out_sz & 3i32 {
            1 => {
                c_0 =
                    D.R[Rans64DecGet(&mut rans0 as *mut RansState,
                                     12i32 as uint32_t) as usize];
                Rans64DecAdvanceSymbol(&mut rans0 as *mut RansState,
                                       &mut ptr as *mut *mut uint8_t as
                                           *mut *mut uint32_t,
                                       &mut syms[c_0 as usize] as
                                           *mut RansDecSymbol,
                                       12i32 as uint32_t);
                *out_buf.offset(out_end as isize) = c_0 as c_char
            }
            2 => {
                c_0 =
                    D.R[Rans64DecGet(&mut rans0 as *mut RansState,
                                     12i32 as uint32_t) as usize];
                Rans64DecAdvanceSymbol(&mut rans0 as *mut RansState,
                                       &mut ptr as *mut *mut uint8_t as
                                           *mut *mut uint32_t,
                                       &mut syms[c_0 as usize] as
                                           *mut RansDecSymbol,
                                       12i32 as uint32_t);
                *out_buf.offset(out_end as isize) = c_0 as c_char;
                c_0 =
                    D.R[Rans64DecGet(&mut rans1 as *mut RansState,
                                     12i32 as uint32_t) as usize];
                Rans64DecAdvanceSymbol(&mut rans1 as *mut RansState,
                                       &mut ptr as *mut *mut uint8_t as
                                           *mut *mut uint32_t,
                                       &mut syms[c_0 as usize] as
                                           *mut RansDecSymbol,
                                       12i32 as uint32_t);
                *out_buf.offset((out_end + 1i32) as isize) =
                    c_0 as c_char
            }
            3 => {
                c_0 =
                    D.R[Rans64DecGet(&mut rans0 as *mut RansState,
                                     12i32 as uint32_t) as usize];
                Rans64DecAdvanceSymbol(&mut rans0 as *mut RansState,
                                       &mut ptr as *mut *mut uint8_t as
                                           *mut *mut uint32_t,
                                       &mut syms[c_0 as usize] as
                                           *mut RansDecSymbol,
                                       12i32 as uint32_t);
                *out_buf.offset(out_end as isize) = c_0 as c_char;
                c_0 =
                    D.R[Rans64DecGet(&mut rans1 as *mut RansState,
                                     12i32 as uint32_t) as usize];
                Rans64DecAdvanceSymbol(&mut rans1 as *mut RansState,
                                       &mut ptr as *mut *mut uint8_t as
                                           *mut *mut uint32_t,
                                       &mut syms[c_0 as usize] as
                                           *mut RansDecSymbol,
                                       12i32 as uint32_t);
                *out_buf.offset((out_end + 1i32) as isize) =
                    c_0 as c_char;
                c_0 =
                    D.R[Rans64DecGet(&mut rans2 as *mut RansState,
                                     12i32 as uint32_t) as usize];
                Rans64DecAdvanceSymbol(&mut rans2 as *mut RansState,
                                       &mut ptr as *mut *mut uint8_t as
                                           *mut *mut uint32_t,
                                       &mut syms[c_0 as usize] as
                                           *mut RansDecSymbol,
                                       12i32 as uint32_t);
                *out_buf.offset((out_end + 2i32) as isize) =
                    c_0 as c_char
            }
            0 | _ => { }
        }
        *out_size = out_sz as c_uint;
        return out_buf as *mut c_uchar
    };
}
unsafe extern "C" fn hist1_4(mut in_0: *mut c_uchar,
                             mut in_size: c_uint,
                             mut F0: *mut [c_int; 256],
                             mut T0: *mut c_int) -> () {
    let mut T1: [c_int; 264] =
        [0i32, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
    let mut T2: [c_int; 264] =
        [0i32, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
    let mut T3: [c_int; 264] =
        [0i32, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
    let mut idiv4: c_uint = in_size.wrapping_div(4i32 as c_uint);
    let mut i: c_int = 0;
    let mut c0: c_uchar = 0;
    let mut c1: c_uchar = 0;
    let mut c2: c_uchar = 0;
    let mut c3: c_uchar = 0;
    let mut in0: *mut c_uchar = in_0.offset(0isize);
    let mut in1: *mut c_uchar = in_0.offset(idiv4 as isize);
    let mut in2: *mut c_uchar =
        in_0.offset(idiv4.wrapping_mul(2i32 as c_uint) as isize);
    let mut in3: *mut c_uchar =
        in_0.offset(idiv4.wrapping_mul(3i32 as c_uint) as isize);
    let mut last_0: c_uchar = 0i32 as c_uchar;
    let mut last_1: c_uchar = *in1.offset(-1i32 as isize);
    let mut last_2: c_uchar = *in2.offset(-1i32 as isize);
    let mut last_3: c_uchar = *in3.offset(-1i32 as isize);
    let mut in0_end: *mut c_uchar = in1;
    while in0 < in0_end {
        let fresh36 = in0;
        in0 = in0.offset(1);
        c0 = *fresh36;
        let ref mut fresh37 = (*F0.offset(last_0 as isize))[c0 as usize];
        *fresh37 += 1;
        let ref mut fresh38 = *T0.offset(last_0 as isize);
        *fresh38 += 1;
        last_0 = c0;
        let fresh39 = in1;
        in1 = in1.offset(1);
        c1 = *fresh39;
        let ref mut fresh40 = (*F0.offset(last_1 as isize))[c1 as usize];
        *fresh40 += 1;
        T1[last_1 as usize] += 1;
        last_1 = c1;
        let fresh41 = in2;
        in2 = in2.offset(1);
        c2 = *fresh41;
        let ref mut fresh42 = (*F0.offset(last_2 as isize))[c2 as usize];
        *fresh42 += 1;
        T2[last_2 as usize] += 1;
        last_2 = c2;
        let fresh43 = in3;
        in3 = in3.offset(1);
        c3 = *fresh43;
        let ref mut fresh44 = (*F0.offset(last_3 as isize))[c3 as usize];
        *fresh44 += 1;
        T3[last_3 as usize] += 1;
        last_3 = c3
    }
    while in3 < in_0.offset(in_size as isize) {
        let fresh45 = in3;
        in3 = in3.offset(1);
        c3 = *fresh45;
        let ref mut fresh46 = (*F0.offset(last_3 as isize))[c3 as usize];
        *fresh46 += 1;
        T3[last_3 as usize] += 1;
        last_3 = c3
    }
    i = 0i32;
    while i < 256i32 {
        *T0.offset(i as isize) +=
            T1[i as usize] + T2[i as usize] + T3[i as usize];
        i += 1
    };
}
#[no_mangle]
pub unsafe extern "C" fn rans_compress_O1(mut in_0: *mut c_uchar,
                                          mut in_size: c_uint,
                                          mut out_size: *mut c_uint)
 -> *mut c_uchar {
    let mut out_buf: *mut c_uchar =
        malloc((1.05f64 * in_size as c_double +
                    (257i32 * 257i32 * 3i32) as c_double +
                    4i32 as c_double) as c_ulong) as
            *mut c_uchar;
    let mut cp: *mut c_uchar = out_buf;
    let mut out_end: *mut c_uchar = 0 as *mut c_uchar;
    let mut tab_size: c_uint = 0;
    let mut rle_i: c_uint = 0;
    let mut rle_j: c_uint = 0;
    let mut syms: [[RansEncSymbol; 256]; 256] =
        [[RansEncSymbol{rcp_freq: 0,
                        freq: 0,
                        bias: 0,
                        cmpl_freq: 0,
                        rcp_shift: 0,}; 256]; 256];
    if out_buf.is_null() {
        return 0 as *mut c_uchar
    } else {
        out_end =
            out_buf.offset((1.05f64 * in_size as c_double) as
                               c_int as
                               isize).offset((257i32 * 257i32 * 3i32) as
                                                 isize).offset(4isize);
        cp = out_buf.offset(4isize);
        let mut F: [[c_int; 256]; 256] =
            [[0i32, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
              0, 0, 0, 0, 0], [0; 256], [0; 256], [0; 256], [0; 256],
             [0; 256], [0; 256], [0; 256], [0; 256], [0; 256], [0; 256],
             [0; 256], [0; 256], [0; 256], [0; 256], [0; 256], [0; 256],
             [0; 256], [0; 256], [0; 256], [0; 256], [0; 256], [0; 256],
             [0; 256], [0; 256], [0; 256], [0; 256], [0; 256], [0; 256],
             [0; 256], [0; 256], [0; 256], [0; 256], [0; 256], [0; 256],
             [0; 256], [0; 256], [0; 256], [0; 256], [0; 256], [0; 256],
             [0; 256], [0; 256], [0; 256], [0; 256], [0; 256], [0; 256],
             [0; 256], [0; 256], [0; 256], [0; 256], [0; 256], [0; 256],
             [0; 256], [0; 256], [0; 256], [0; 256], [0; 256], [0; 256],
             [0; 256], [0; 256], [0; 256], [0; 256], [0; 256], [0; 256],
             [0; 256], [0; 256], [0; 256], [0; 256], [0; 256], [0; 256],
             [0; 256], [0; 256], [0; 256], [0; 256], [0; 256], [0; 256],
             [0; 256], [0; 256], [0; 256], [0; 256], [0; 256], [0; 256],
             [0; 256], [0; 256], [0; 256], [0; 256], [0; 256], [0; 256],
             [0; 256], [0; 256], [0; 256], [0; 256], [0; 256], [0; 256],
             [0; 256], [0; 256], [0; 256], [0; 256], [0; 256], [0; 256],
             [0; 256], [0; 256], [0; 256], [0; 256], [0; 256], [0; 256],
             [0; 256], [0; 256], [0; 256], [0; 256], [0; 256], [0; 256],
             [0; 256], [0; 256], [0; 256], [0; 256], [0; 256], [0; 256],
             [0; 256], [0; 256], [0; 256], [0; 256], [0; 256], [0; 256],
             [0; 256], [0; 256], [0; 256], [0; 256], [0; 256], [0; 256],
             [0; 256], [0; 256], [0; 256], [0; 256], [0; 256], [0; 256],
             [0; 256], [0; 256], [0; 256], [0; 256], [0; 256], [0; 256],
             [0; 256], [0; 256], [0; 256], [0; 256], [0; 256], [0; 256],
             [0; 256], [0; 256], [0; 256], [0; 256], [0; 256], [0; 256],
             [0; 256], [0; 256], [0; 256], [0; 256], [0; 256], [0; 256],
             [0; 256], [0; 256], [0; 256], [0; 256], [0; 256], [0; 256],
             [0; 256], [0; 256], [0; 256], [0; 256], [0; 256], [0; 256],
             [0; 256], [0; 256], [0; 256], [0; 256], [0; 256], [0; 256],
             [0; 256], [0; 256], [0; 256], [0; 256], [0; 256], [0; 256],
             [0; 256], [0; 256], [0; 256], [0; 256], [0; 256], [0; 256],
             [0; 256], [0; 256], [0; 256], [0; 256], [0; 256], [0; 256],
             [0; 256], [0; 256], [0; 256], [0; 256], [0; 256], [0; 256],
             [0; 256], [0; 256], [0; 256], [0; 256], [0; 256], [0; 256],
             [0; 256], [0; 256], [0; 256], [0; 256], [0; 256], [0; 256],
             [0; 256], [0; 256], [0; 256], [0; 256], [0; 256], [0; 256],
             [0; 256], [0; 256], [0; 256], [0; 256], [0; 256], [0; 256],
             [0; 256], [0; 256], [0; 256], [0; 256], [0; 256], [0; 256],
             [0; 256], [0; 256], [0; 256], [0; 256], [0; 256], [0; 256],
             [0; 256], [0; 256], [0; 256], [0; 256], [0; 256], [0; 256],
             [0; 256], [0; 256], [0; 256], [0; 256], [0; 256], [0; 256],
             [0; 256], [0; 256], [0; 256], [0; 256], [0; 256]];
        let mut T: [c_int; 264] =
            [0i32, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
             0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
             0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
             0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
             0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
             0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
             0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
             0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
             0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
             0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
             0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
             0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
             0];
        let mut i: c_int = 0;
        let mut j: c_int = 0;
        hist1_4(in_0, in_size, F.as_mut_ptr(), T.as_mut_ptr());
        F[0usize][*in_0.offset((1i32 as
                                    c_uint).wrapping_mul(in_size >>
                                                                   2i32) as
                                   isize) as usize] += 1;
        F[0usize][*in_0.offset((2i32 as
                                    c_uint).wrapping_mul(in_size >>
                                                                   2i32) as
                                   isize) as usize] += 1;
        F[0usize][*in_0.offset((3i32 as
                                    c_uint).wrapping_mul(in_size >>
                                                                   2i32) as
                                   isize) as usize] += 1;
        T[0usize] += 3i32;
        i = 0i32;
        rle_i = i as c_uint;
        while i < 256i32 {
            let mut t2: c_int = 0;
            let mut m: c_int = 0;
            let mut M: c_int = 0;
            let mut x: c_uint = 0;
            if !(T[i as usize] == 0i32) {
                let mut p: c_double =
                    (1i32 << 12i32) as c_double /
                        T[i as usize] as c_double;
                j = 0i32;
                M = j;
                m = M;
                t2 = m;
                while j < 256i32 {
                    if !(0 == F[i as usize][j as usize]) {
                        if m < F[i as usize][j as usize] {
                            m = F[i as usize][j as usize];
                            M = j
                        }
                        F[i as usize][j as usize] =
                            (F[i as usize][j as usize] as c_double * p)
                                as c_int;
                        if F[i as usize][j as usize] == 0i32 {
                            F[i as usize][j as usize] = 1i32
                        }
                        t2 += F[i as usize][j as usize]
                    }
                    j += 1
                }
                t2 += 1;
                if t2 < 1i32 << 12i32 {
                    F[i as usize][M as usize] += (1i32 << 12i32) - t2
                } else { F[i as usize][M as usize] -= t2 - (1i32 << 12i32) }
                if 0 != rle_i {
                    rle_i = rle_i.wrapping_sub(1)
                } else {
                    let fresh47 = cp;
                    cp = cp.offset(1);
                    *fresh47 = i as c_uchar;
                    if 0 != i && 0 != T[(i - 1i32) as usize] {
                        rle_i = (i + 1i32) as c_uint;
                        while rle_i < 256i32 as c_uint &&
                                  0 != T[rle_i as usize] {
                            rle_i = rle_i.wrapping_add(1)
                        }
                        rle_i =
                            rle_i.wrapping_sub((i + 1i32) as c_uint);
                        let fresh48 = cp;
                        cp = cp.offset(1);
                        *fresh48 = rle_i as c_uchar
                    }
                }
                let mut F_i_: *mut c_int = F[i as usize].as_mut_ptr();
                x = 0i32 as c_uint;
                rle_j = 0i32 as c_uint;
                j = 0i32;
                while j < 256i32 {
                    if 0 != *F_i_.offset(j as isize) {
                        if 0 != rle_j {
                            rle_j = rle_j.wrapping_sub(1)
                        } else {
                            let fresh49 = cp;
                            cp = cp.offset(1);
                            *fresh49 = j as c_uchar;
                            if 0 == rle_j && 0 != j &&
                                   0 != *F_i_.offset((j - 1i32) as isize) {
                                rle_j = (j + 1i32) as c_uint;
                                while rle_j < 256i32 as c_uint &&
                                          0 != *F_i_.offset(rle_j as isize) {
                                    rle_j = rle_j.wrapping_add(1)
                                }
                                rle_j =
                                    rle_j.wrapping_sub((j + 1i32) as
                                                           c_uint);
                                let fresh50 = cp;
                                cp = cp.offset(1);
                                *fresh50 = rle_j as c_uchar
                            }
                        }
                        if *F_i_.offset(j as isize) < 128i32 {
                            let fresh51 = cp;
                            cp = cp.offset(1);
                            *fresh51 =
                                *F_i_.offset(j as isize) as c_uchar
                        } else {
                            let fresh52 = cp;
                            cp = cp.offset(1);
                            *fresh52 =
                                (128i32 | *F_i_.offset(j as isize) >> 8i32) as
                                    c_uchar;
                            let fresh53 = cp;
                            cp = cp.offset(1);
                            *fresh53 =
                                (*F_i_.offset(j as isize) & 255i32) as
                                    c_uchar
                        }
                        Rans64EncSymbolInit(&mut syms[i as usize][j as usize]
                                                as *mut RansEncSymbol, x,
                                            *F_i_.offset(j as isize) as
                                                uint32_t, 12i32 as uint32_t);
                        x =
                            x.wrapping_add(*F_i_.offset(j as isize) as
                                               c_uint)
                    }
                    j += 1
                }
                let fresh54 = cp;
                cp = cp.offset(1);
                *fresh54 = 0i32 as c_uchar
            }
            i += 1
        }
        let fresh55 = cp;
        cp = cp.offset(1);
        *fresh55 = 0i32 as c_uchar;
        tab_size =
            out_buf.offset_to(cp).expect("bad offset_to") as c_long as
                c_uint;
        let mut rans0: RansState = 0;
        let mut rans1: RansState = 0;
        let mut rans2: RansState = 0;
        let mut rans3: RansState = 0;
        Rans64EncInit(&mut rans0 as *mut RansState);
        Rans64EncInit(&mut rans1 as *mut RansState);
        Rans64EncInit(&mut rans2 as *mut RansState);
        Rans64EncInit(&mut rans3 as *mut RansState);
        let mut ptr: *mut uint8_t = out_end;
        let mut isz4: c_int = (in_size >> 2i32) as c_int;
        let mut i0: c_int = 1i32 * isz4 - 2i32;
        let mut i1: c_int = 2i32 * isz4 - 2i32;
        let mut i2: c_int = 3i32 * isz4 - 2i32;
        let mut i3: c_int = 4i32 * isz4 - 2i32;
        let mut l0: c_uchar = *in_0.offset((i0 + 1i32) as isize);
        let mut l1: c_uchar = *in_0.offset((i1 + 1i32) as isize);
        let mut l2: c_uchar = *in_0.offset((i2 + 1i32) as isize);
        let mut l3: c_uchar = *in_0.offset((i3 + 1i32) as isize);
        l3 =
            *in_0.offset(in_size.wrapping_sub(1i32 as c_uint) as isize);
        i3 = in_size.wrapping_sub(2i32 as c_uint) as c_int;
        while i3 > 4i32 * isz4 - 2i32 {
            let mut c3: c_uchar = *in_0.offset(i3 as isize);
            Rans64EncPutSymbol(&mut rans3 as *mut RansState,
                               &mut ptr as *mut *mut uint8_t as
                                   *mut *mut uint32_t,
                               &mut syms[c3 as usize][l3 as usize] as
                                   *mut RansEncSymbol, 12i32 as uint32_t);
            l3 = c3;
            i3 -= 1
        }
        while i0 >= 0i32 {
            let mut c0: c_uchar = 0;
            let mut c1: c_uchar = 0;
            let mut c2: c_uchar = 0;
            let mut c3_0: c_uchar = 0;
            c3_0 = *in_0.offset(i3 as isize);
            let mut s3: *mut RansEncSymbol =
                &mut syms[c3_0 as usize][l3 as usize] as *mut RansEncSymbol;
            c2 = *in_0.offset(i2 as isize);
            let mut s2: *mut RansEncSymbol =
                &mut syms[c2 as usize][l2 as usize] as *mut RansEncSymbol;
            c1 = *in_0.offset(i1 as isize);
            let mut s1: *mut RansEncSymbol =
                &mut syms[c1 as usize][l1 as usize] as *mut RansEncSymbol;
            c0 = *in_0.offset(i0 as isize);
            let mut s0: *mut RansEncSymbol =
                &mut syms[c0 as usize][l0 as usize] as *mut RansEncSymbol;
            Rans64EncPutSymbol(&mut rans3 as *mut RansState,
                               &mut ptr as *mut *mut uint8_t as
                                   *mut *mut uint32_t, s3, 12i32 as uint32_t);
            Rans64EncPutSymbol(&mut rans2 as *mut RansState,
                               &mut ptr as *mut *mut uint8_t as
                                   *mut *mut uint32_t, s2, 12i32 as uint32_t);
            Rans64EncPutSymbol(&mut rans1 as *mut RansState,
                               &mut ptr as *mut *mut uint8_t as
                                   *mut *mut uint32_t, s1, 12i32 as uint32_t);
            Rans64EncPutSymbol(&mut rans0 as *mut RansState,
                               &mut ptr as *mut *mut uint8_t as
                                   *mut *mut uint32_t, s0, 12i32 as uint32_t);
            l0 = c0;
            l1 = c1;
            l2 = c2;
            l3 = c3_0;
            i0 -= 1;
            i1 -= 1;
            i2 -= 1;
            i3 -= 1
        }
        Rans64EncPutSymbol(&mut rans3 as *mut RansState,
                           &mut ptr as *mut *mut uint8_t as
                               *mut *mut uint32_t,
                           &mut syms[0usize][l3 as usize] as
                               *mut RansEncSymbol, 12i32 as uint32_t);
        Rans64EncPutSymbol(&mut rans2 as *mut RansState,
                           &mut ptr as *mut *mut uint8_t as
                               *mut *mut uint32_t,
                           &mut syms[0usize][l2 as usize] as
                               *mut RansEncSymbol, 12i32 as uint32_t);
        Rans64EncPutSymbol(&mut rans1 as *mut RansState,
                           &mut ptr as *mut *mut uint8_t as
                               *mut *mut uint32_t,
                           &mut syms[0usize][l1 as usize] as
                               *mut RansEncSymbol, 12i32 as uint32_t);
        Rans64EncPutSymbol(&mut rans0 as *mut RansState,
                           &mut ptr as *mut *mut uint8_t as
                               *mut *mut uint32_t,
                           &mut syms[0usize][l0 as usize] as
                               *mut RansEncSymbol, 12i32 as uint32_t);
        Rans64EncFlush(&mut rans3 as *mut RansState,
                       &mut ptr as *mut *mut uint8_t as *mut *mut uint32_t);
        Rans64EncFlush(&mut rans2 as *mut RansState,
                       &mut ptr as *mut *mut uint8_t as *mut *mut uint32_t);
        Rans64EncFlush(&mut rans1 as *mut RansState,
                       &mut ptr as *mut *mut uint8_t as *mut *mut uint32_t);
        Rans64EncFlush(&mut rans0 as *mut RansState,
                       &mut ptr as *mut *mut uint8_t as *mut *mut uint32_t);
        *out_size =
            (ptr.offset_to(out_end).expect("bad offset_to") as c_long +
                 tab_size as c_long) as c_uint;
        cp = out_buf;
        let fresh56 = cp;
        cp = cp.offset(1);
        *fresh56 =
            (in_size >> 0i32 & 255i32 as c_uint) as c_uchar;
        let fresh57 = cp;
        cp = cp.offset(1);
        *fresh57 =
            (in_size >> 8i32 & 255i32 as c_uint) as c_uchar;
        let fresh58 = cp;
        cp = cp.offset(1);
        *fresh58 =
            (in_size >> 16i32 & 255i32 as c_uint) as c_uchar;
        let fresh59 = cp;
        cp = cp.offset(1);
        *fresh59 =
            (in_size >> 24i32 & 255i32 as c_uint) as c_uchar;
        memmove(out_buf.offset(tab_size as isize) as *mut c_void,
                ptr as *const c_void,
                ptr.offset_to(out_end).expect("bad offset_to") as c_long
                    as c_ulong);
        return out_buf
    };
}
#[no_mangle]
pub unsafe extern "C" fn rans_uncompress_O1(mut in_0: *mut c_uchar,
                                            mut in_size: c_uint,
                                            mut out_size: *mut c_uint)
 -> *mut c_uchar {
    let mut cp: *mut c_uchar = in_0.offset(4isize);
    let mut i: c_int = 0;
    let mut j: c_int = -999i32;
    let mut x: c_int = 0;
    let mut out_sz: c_int = 0;
    let mut rle_i: c_int = 0;
    let mut rle_j: c_int = 0;
    let mut out_buf: *mut c_char = 0 as *mut c_char;
    let mut D: [ari_decoder; 256] = [ari_decoder{R: [0; 4096],}; 256];
    let mut syms: [[RansDecSymbol; 256]; 256] =
        [[RansDecSymbol{start: 0, freq: 0,}; 256]; 256];
    out_sz =
        (*in_0.offset(0isize) as c_int) << 0i32 |
            (*in_0.offset(1isize) as c_int) << 8i32 |
            (*in_0.offset(2isize) as c_int) << 16i32 |
            (*in_0.offset(3isize) as c_int) << 24i32;
    out_buf = malloc(out_sz as c_ulong) as *mut c_char;
    if out_buf.is_null() {
        return 0 as *mut c_uchar
    } else {
        rle_i = 0i32;
        let fresh60 = cp;
        cp = cp.offset(1);
        i = *fresh60 as c_int;
        loop  {
            x = 0i32;
            rle_j = x;
            let fresh61 = cp;
            cp = cp.offset(1);
            j = *fresh61 as c_int;
            loop  {
                let mut F: c_int = 0;
                let mut C: c_int = 0;
                let fresh62 = cp;
                cp = cp.offset(1);
                F = *fresh62 as c_int;
                if F >= 128i32 {
                    F &= !128i32;
                    let fresh63 = cp;
                    cp = cp.offset(1);
                    F = (F & 127i32) << 8i32 | *fresh63 as c_int
                }
                C = x;
                if 0 == F { F = 1i32 << 12i32 }
                Rans64DecSymbolInit(&mut syms[i as usize][j as usize] as
                                        *mut RansDecSymbol, C as uint32_t,
                                    F as uint32_t);
                memset(&mut D[i as usize].R[x as usize] as *mut c_uchar
                           as *mut c_void, j, F as c_ulong);
                x += F;
                if 0 == rle_j && j + 1i32 == *cp as c_int {
                    let fresh64 = cp;
                    cp = cp.offset(1);
                    j = *fresh64 as c_int;
                    let fresh65 = cp;
                    cp = cp.offset(1);
                    rle_j = *fresh65 as c_int
                } else if 0 != rle_j {
                    rle_j -= 1;
                    j += 1
                } else {
                    let fresh66 = cp;
                    cp = cp.offset(1);
                    j = *fresh66 as c_int
                }
                if !(0 != j) { break ; }
            }
            if 0 == rle_i && i + 1i32 == *cp as c_int {
                let fresh67 = cp;
                cp = cp.offset(1);
                i = *fresh67 as c_int;
                let fresh68 = cp;
                cp = cp.offset(1);
                rle_i = *fresh68 as c_int
            } else if 0 != rle_i {
                rle_i -= 1;
                i += 1
            } else {
                let fresh69 = cp;
                cp = cp.offset(1);
                i = *fresh69 as c_int
            }
            if !(0 != i) { break ; }
        }
        let mut rans0: RansState = 0;
        let mut rans1: RansState = 0;
        let mut rans2: RansState = 0;
        let mut rans3: RansState = 0;
        let mut ptr: *mut uint8_t = cp;
        Rans64DecInit(&mut rans0 as *mut RansState,
                      &mut ptr as *mut *mut uint8_t as *mut *mut uint32_t);
        Rans64DecInit(&mut rans1 as *mut RansState,
                      &mut ptr as *mut *mut uint8_t as *mut *mut uint32_t);
        Rans64DecInit(&mut rans2 as *mut RansState,
                      &mut ptr as *mut *mut uint8_t as *mut *mut uint32_t);
        Rans64DecInit(&mut rans3 as *mut RansState,
                      &mut ptr as *mut *mut uint8_t as *mut *mut uint32_t);
        let mut isz4: c_int = out_sz >> 2i32;
        let mut l0: c_int = 0i32;
        let mut l1: c_int = 0i32;
        let mut l2: c_int = 0i32;
        let mut l3: c_int = 0i32;
        let mut i4: [c_int; 4] =
            [0i32 * isz4, 1i32 * isz4, 2i32 * isz4, 3i32 * isz4];
        let mut R: [RansState; 4] = [0; 4];
        R[0usize] = rans0;
        R[1usize] = rans1;
        R[2usize] = rans2;
        R[3usize] = rans3;
        while i4[0usize] < isz4 {
            let mut m: [uint32_t; 4] =
                [(R[0usize] &
                      (1u32 << 12i32).wrapping_sub(1i32 as c_uint) as
                          c_ulong) as uint32_t,
                 (R[1usize] &
                      (1u32 << 12i32).wrapping_sub(1i32 as c_uint) as
                          c_ulong) as uint32_t,
                 (R[2usize] &
                      (1u32 << 12i32).wrapping_sub(1i32 as c_uint) as
                          c_ulong) as uint32_t,
                 (R[3usize] &
                      (1u32 << 12i32).wrapping_sub(1i32 as c_uint) as
                          c_ulong) as uint32_t];
            let mut c: [uint8_t; 4] =
                [D[l0 as usize].R[m[0usize] as usize],
                 D[l1 as usize].R[m[1usize] as usize],
                 D[l2 as usize].R[m[2usize] as usize],
                 D[l3 as usize].R[m[3usize] as usize]];
            *out_buf.offset(i4[0usize] as isize) = c[0usize] as c_char;
            *out_buf.offset(i4[1usize] as isize) = c[1usize] as c_char;
            *out_buf.offset(i4[2usize] as isize) = c[2usize] as c_char;
            *out_buf.offset(i4[3usize] as isize) = c[3usize] as c_char;
            R[0usize] =
                (syms[l0 as usize][c[0usize] as usize].freq as
                     c_ulong).wrapping_mul(R[0usize] >> 12i32);
            R[1usize] =
                (syms[l1 as usize][c[1usize] as usize].freq as
                     c_ulong).wrapping_mul(R[1usize] >> 12i32);
            R[2usize] =
                (syms[l2 as usize][c[2usize] as usize].freq as
                     c_ulong).wrapping_mul(R[2usize] >> 12i32);
            R[3usize] =
                (syms[l3 as usize][c[3usize] as usize].freq as
                     c_ulong).wrapping_mul(R[3usize] >> 12i32);
            R[0usize] =
                (R[0usize] as
                     c_ulong).wrapping_add(m[0usize].wrapping_sub(syms[l0
                                                                                 as
                                                                                 usize][c[0usize]
                                                                                            as
                                                                                            usize].start)
                                                     as c_ulong) as
                    RansState as RansState;
            R[1usize] =
                (R[1usize] as
                     c_ulong).wrapping_add(m[1usize].wrapping_sub(syms[l1
                                                                                 as
                                                                                 usize][c[1usize]
                                                                                            as
                                                                                            usize].start)
                                                     as c_ulong) as
                    RansState as RansState;
            R[2usize] =
                (R[2usize] as
                     c_ulong).wrapping_add(m[2usize].wrapping_sub(syms[l2
                                                                                 as
                                                                                 usize][c[2usize]
                                                                                            as
                                                                                            usize].start)
                                                     as c_ulong) as
                    RansState as RansState;
            R[3usize] =
                (R[3usize] as
                     c_ulong).wrapping_add(m[3usize].wrapping_sub(syms[l3
                                                                                 as
                                                                                 usize][c[3usize]
                                                                                            as
                                                                                            usize].start)
                                                     as c_ulong) as
                    RansState as RansState;
            Rans64DecRenorm(&mut R[0usize] as *mut RansState,
                            &mut ptr as *mut *mut uint8_t as
                                *mut *mut uint32_t);
            Rans64DecRenorm(&mut R[1usize] as *mut RansState,
                            &mut ptr as *mut *mut uint8_t as
                                *mut *mut uint32_t);
            Rans64DecRenorm(&mut R[2usize] as *mut RansState,
                            &mut ptr as *mut *mut uint8_t as
                                *mut *mut uint32_t);
            Rans64DecRenorm(&mut R[3usize] as *mut RansState,
                            &mut ptr as *mut *mut uint8_t as
                                *mut *mut uint32_t);
            l0 = c[0usize] as c_int;
            l1 = c[1usize] as c_int;
            l2 = c[2usize] as c_int;
            l3 = c[3usize] as c_int;
            i4[0usize] += 1;
            i4[1usize] += 1;
            i4[2usize] += 1;
            i4[3usize] += 1
        }
        rans0 = R[0usize];
        rans1 = R[1usize];
        rans2 = R[2usize];
        rans3 = R[3usize];
        while i4[3usize] < out_sz {
            let mut c3: c_uchar =
                D[l3 as
                      usize].R[Rans64DecGet(&mut rans3 as *mut RansState,
                                            12i32 as uint32_t) as usize];
            *out_buf.offset(i4[3usize] as isize) = c3 as c_char;
            Rans64DecAdvanceSymbol(&mut rans3 as *mut RansState,
                                   &mut ptr as *mut *mut uint8_t as
                                       *mut *mut uint32_t,
                                   &mut syms[l3 as usize][c3 as usize] as
                                       *mut RansDecSymbol, 12i32 as uint32_t);
            l3 = c3 as c_int;
            i4[3usize] += 1
        }
        *out_size = out_sz as c_uint;
        return out_buf as *mut c_uchar
    };
}
#[no_mangle]
pub unsafe extern "C" fn rans_uncompress_O1b(mut in_0: *mut c_uchar,
                                             mut in_size: c_uint,
                                             mut out_size: *mut c_uint)
 -> *mut c_uchar {
    let mut cp: *mut c_uchar = in_0.offset(4isize);
    let mut i: c_int = 0;
    let mut j: c_int = -999i32;
    let mut x: c_int = 0;
    let mut out_sz: c_int = 0;
    let mut rle_i: c_int = 0;
    let mut rle_j: c_int = 0;
    let mut out_buf: *mut c_char = 0 as *mut c_char;
    let mut D: [ari_decoder; 256] = [ari_decoder{R: [0; 4096],}; 256];
    let mut syms: [[RansDecSymbol; 256]; 256] =
        [[RansDecSymbol{start: 0, freq: 0,}; 256]; 256];
    out_sz =
        (*in_0.offset(0isize) as c_int) << 0i32 |
            (*in_0.offset(1isize) as c_int) << 8i32 |
            (*in_0.offset(2isize) as c_int) << 16i32 |
            (*in_0.offset(3isize) as c_int) << 24i32;
    out_buf = malloc(out_sz as c_ulong) as *mut c_char;
    if out_buf.is_null() {
        return 0 as *mut c_uchar
    } else {
        rle_i = 0i32;
        let fresh70 = cp;
        cp = cp.offset(1);
        i = *fresh70 as c_int;
        loop  {
            x = 0i32;
            rle_j = x;
            let fresh71 = cp;
            cp = cp.offset(1);
            j = *fresh71 as c_int;
            loop  {
                let mut F: c_int = 0;
                let mut C: c_int = 0;
                let fresh72 = cp;
                cp = cp.offset(1);
                F = *fresh72 as c_int;
                if F >= 128i32 {
                    F &= !128i32;
                    let fresh73 = cp;
                    cp = cp.offset(1);
                    F = (F & 127i32) << 8i32 | *fresh73 as c_int
                }
                C = x;
                if 0 == F { F = 1i32 << 12i32 }
                Rans64DecSymbolInit(&mut syms[i as usize][j as usize] as
                                        *mut RansDecSymbol, C as uint32_t,
                                    F as uint32_t);
                memset(&mut D[i as usize].R[x as usize] as *mut c_uchar
                           as *mut c_void, j, F as c_ulong);
                x += F;
                if 0 == rle_j && j + 1i32 == *cp as c_int {
                    let fresh74 = cp;
                    cp = cp.offset(1);
                    j = *fresh74 as c_int;
                    let fresh75 = cp;
                    cp = cp.offset(1);
                    rle_j = *fresh75 as c_int
                } else if 0 != rle_j {
                    rle_j -= 1;
                    j += 1
                } else {
                    let fresh76 = cp;
                    cp = cp.offset(1);
                    j = *fresh76 as c_int
                }
                if !(0 != j) { break ; }
            }
            if 0 == rle_i && i + 1i32 == *cp as c_int {
                let fresh77 = cp;
                cp = cp.offset(1);
                i = *fresh77 as c_int;
                let fresh78 = cp;
                cp = cp.offset(1);
                rle_i = *fresh78 as c_int
            } else if 0 != rle_i {
                rle_i -= 1;
                i += 1
            } else {
                let fresh79 = cp;
                cp = cp.offset(1);
                i = *fresh79 as c_int
            }
            if !(0 != i) { break ; }
        }
        let mut rans0: RansState = 0;
        let mut rans1: RansState = 0;
        let mut rans2: RansState = 0;
        let mut rans3: RansState = 0;
        let mut ptr: *mut uint8_t = cp;
        Rans64DecInit(&mut rans0 as *mut RansState,
                      &mut ptr as *mut *mut uint8_t as *mut *mut uint32_t);
        Rans64DecInit(&mut rans1 as *mut RansState,
                      &mut ptr as *mut *mut uint8_t as *mut *mut uint32_t);
        Rans64DecInit(&mut rans2 as *mut RansState,
                      &mut ptr as *mut *mut uint8_t as *mut *mut uint32_t);
        Rans64DecInit(&mut rans3 as *mut RansState,
                      &mut ptr as *mut *mut uint8_t as *mut *mut uint32_t);
        let mut isz4: c_int = out_sz >> 2i32;
        let mut l0: c_int = 0i32;
        let mut l1: c_int = 0i32;
        let mut l2: c_int = 0i32;
        let mut l3: c_int = 0i32;
        let mut i4: [c_int; 4] =
            [0i32 * isz4, 1i32 * isz4, 2i32 * isz4, 3i32 * isz4];
        let mut R: [RansState; 4] = [0; 4];
        R[0usize] = rans0;
        R[1usize] = rans1;
        R[2usize] = rans2;
        R[3usize] = rans3;
        while i4[0usize] < isz4 {
            let mut m: [uint32_t; 4] =
                [(R[0usize] &
                      (1u32 << 12i32).wrapping_sub(1i32 as c_uint) as
                          c_ulong) as uint32_t,
                 (R[1usize] &
                      (1u32 << 12i32).wrapping_sub(1i32 as c_uint) as
                          c_ulong) as uint32_t,
                 (R[2usize] &
                      (1u32 << 12i32).wrapping_sub(1i32 as c_uint) as
                          c_ulong) as uint32_t,
                 (R[3usize] &
                      (1u32 << 12i32).wrapping_sub(1i32 as c_uint) as
                          c_ulong) as uint32_t];
            let mut c: [uint8_t; 4] =
                [D[l0 as usize].R[m[0usize] as usize],
                 D[l1 as usize].R[m[1usize] as usize],
                 D[l2 as usize].R[m[2usize] as usize],
                 D[l3 as usize].R[m[3usize] as usize]];
            *out_buf.offset(i4[0usize] as isize) = c[0usize] as c_char;
            *out_buf.offset(i4[1usize] as isize) = c[1usize] as c_char;
            *out_buf.offset(i4[2usize] as isize) = c[2usize] as c_char;
            *out_buf.offset(i4[3usize] as isize) = c[3usize] as c_char;
            R[0usize] =
                (syms[l0 as usize][c[0usize] as usize].freq as
                     c_ulong).wrapping_mul(R[0usize] >> 12i32);
            R[0usize] =
                (R[0usize] as
                     c_ulong).wrapping_add(m[0usize].wrapping_sub(syms[l0
                                                                                 as
                                                                                 usize][c[0usize]
                                                                                            as
                                                                                            usize].start)
                                                     as c_ulong) as
                    RansState as RansState;
            R[1usize] =
                (syms[l1 as usize][c[1usize] as usize].freq as
                     c_ulong).wrapping_mul(R[1usize] >> 12i32);
            R[1usize] =
                (R[1usize] as
                     c_ulong).wrapping_add(m[1usize].wrapping_sub(syms[l1
                                                                                 as
                                                                                 usize][c[1usize]
                                                                                            as
                                                                                            usize].start)
                                                     as c_ulong) as
                    RansState as RansState;
            R[2usize] =
                (syms[l2 as usize][c[2usize] as usize].freq as
                     c_ulong).wrapping_mul(R[2usize] >> 12i32);
            R[2usize] =
                (R[2usize] as
                     c_ulong).wrapping_add(m[2usize].wrapping_sub(syms[l2
                                                                                 as
                                                                                 usize][c[2usize]
                                                                                            as
                                                                                            usize].start)
                                                     as c_ulong) as
                    RansState as RansState;
            R[3usize] =
                (syms[l3 as usize][c[3usize] as usize].freq as
                     c_ulong).wrapping_mul(R[3usize] >> 12i32);
            R[3usize] =
                (R[3usize] as
                     c_ulong).wrapping_add(m[3usize].wrapping_sub(syms[l3
                                                                                 as
                                                                                 usize][c[3usize]
                                                                                            as
                                                                                            usize].start)
                                                     as c_ulong) as
                    RansState as RansState;
            Rans64DecRenorm(&mut R[0usize] as *mut RansState,
                            &mut ptr as *mut *mut uint8_t as
                                *mut *mut uint32_t);
            Rans64DecRenorm(&mut R[1usize] as *mut RansState,
                            &mut ptr as *mut *mut uint8_t as
                                *mut *mut uint32_t);
            Rans64DecRenorm(&mut R[2usize] as *mut RansState,
                            &mut ptr as *mut *mut uint8_t as
                                *mut *mut uint32_t);
            Rans64DecRenorm(&mut R[3usize] as *mut RansState,
                            &mut ptr as *mut *mut uint8_t as
                                *mut *mut uint32_t);
            l0 = c[0usize] as c_int;
            l1 = c[1usize] as c_int;
            l2 = c[2usize] as c_int;
            l3 = c[3usize] as c_int;
            i4[0usize] += 1;
            i4[1usize] += 1;
            i4[2usize] += 1;
            i4[3usize] += 1
        }
        rans0 = R[0usize];
        rans1 = R[1usize];
        rans2 = R[2usize];
        rans3 = R[3usize];
        while i4[3usize] < out_sz {
            let mut c3: c_uchar =
                D[l3 as
                      usize].R[Rans64DecGet(&mut rans3 as *mut RansState,
                                            12i32 as uint32_t) as usize];
            *out_buf.offset(i4[3usize] as isize) = c3 as c_char;
            Rans64DecAdvanceSymbol(&mut rans3 as *mut RansState,
                                   &mut ptr as *mut *mut uint8_t as
                                       *mut *mut uint32_t,
                                   &mut syms[l3 as usize][c3 as usize] as
                                       *mut RansDecSymbol, 12i32 as uint32_t);
            l3 = c3 as c_int;
            i4[3usize] += 1
        }
        *out_size = out_sz as c_uint;
        return out_buf as *mut c_uchar
    };
}
#[no_mangle]
pub unsafe extern "C" fn rans_compress(mut in_0: *mut c_uchar,
                                       mut in_size: c_uint,
                                       mut out_size: *mut c_uint,
                                       mut order: c_int)
 -> *mut c_uchar {
    return if 0 != order {
               rans_compress_O1(in_0, in_size, out_size)
           } else { rans_compress_O0(in_0, in_size, out_size) };
}
#[no_mangle]
pub unsafe extern "C" fn rans_uncompress(mut in_0: *mut c_uchar,
                                         mut in_size: c_uint,
                                         mut out_size: *mut c_uint,
                                         mut order: c_int)
 -> *mut c_uchar {
    return if 0 != order {
               rans_uncompress_O1b(in_0, in_size, out_size)
           } else { rans_uncompress_O0b(in_0, in_size, out_size) };
}

pub fn main() {
    unsafe{main_0()};
}
pub unsafe extern "C" fn main_0() -> c_int {
    let mut opt: c_int = 0;
    let mut order: c_int = 0i32;
    let mut in_buf: [c_uchar; 1299151] = [0; 1299151];
    let mut decode: c_int = 0i32;
    let mut test: c_int = 1i32;
    let mut tv1: timeval = timeval{tv_sec: 0, tv_usec: 0,};
    let mut tv2: timeval = timeval{tv_sec: 0, tv_usec: 0,};
    let mut tv3: timeval = timeval{tv_sec: 0, tv_usec: 0,};
    let mut bytes: size_t = 0i32 as size_t;
    extern "C" {
        #[link_name = "optarg"]
        static mut optarg_0: *mut c_char;
    }
    extern "C" {
        #[link_name = "optind"]
        static mut optind_0: c_int;
    }
    order = if 0 != order { 1i32 } else { 0i32 };
    gettimeofday(&mut tv1 as *mut timeval, 0 as *mut timezone);
    if 0 != test {
        let mut len: size_t = 0;
        let mut in_sz: size_t = 0i32 as size_t;
        let mut out_sz: size_t = 0i32 as size_t;
        let mut xlen: size_t = len;
        let mut all_blocks: [blocks; 1] =
            [blocks{blk: 0 as *mut c_uchar, sz: 0,}; 1];
        let mut b: *mut blocks = &mut all_blocks[0usize] as *mut blocks;
        let mut bc: *mut blocks = 0 as *mut blocks;
        let mut bu: *mut blocks = 0 as *mut blocks;
        let mut nb: c_int = 0i32;
        let mut i: c_int = 0;
        len = 1000000i32 as size_t;
        let ref mut fresh80 = (*b.offset(nb as isize)).blk;
        *fresh80 = malloc(len) as *mut c_uchar;
        (*b.offset(nb as isize)).sz = len as uint32_t;
        xlen = len;
        while 0 != xlen {
            xlen = xlen.wrapping_sub(1);
            *(*b.offset(nb as isize)).blk.offset(xlen as isize) =
                (xlen & 63i32 as c_ulong) as c_uchar
        }
        nb += 1;
        in_sz =
            (in_sz as c_ulong).wrapping_add(len) as size_t as size_t;
        let mut trials: c_int = 10i32;
        loop  {
            let fresh81 = trials;
            trials = trials - 1;
            if !(0 != fresh81) { break ; }
            bc =
                malloc((nb as
                            c_ulong).wrapping_mul(::std::mem::size_of::<blocks>()
                                                            as c_ulong))
                    as *mut blocks;
            bu =
                malloc((nb as
                            c_ulong).wrapping_mul(::std::mem::size_of::<blocks>()
                                                            as c_ulong))
                    as *mut blocks;
            gettimeofday(&mut tv1 as *mut timeval, 0 as *mut timezone);
            out_sz = 0i32 as size_t;
            i = 0i32;
            while i < nb {
                let ref mut fresh82 = (*bc.offset(i as isize)).blk;
                *fresh82 =
                    rans_compress((*b.offset(i as isize)).blk,
                                  (*b.offset(i as isize)).sz,
                                  &mut (*bc.offset(i as isize)).sz as
                                      *mut uint32_t, order);
                out_sz =
                    (out_sz as
                         c_ulong).wrapping_add((5i32 as
                                                          c_uint).wrapping_add((*bc.offset(i
                                                                                                     as
                                                                                                     isize)).sz)
                                                         as c_ulong) as
                        size_t as size_t;
                i += 1
            }
            gettimeofday(&mut tv2 as *mut timeval, 0 as *mut timezone);
            i = 0i32;
            while i < nb {
                let ref mut fresh83 = (*bu.offset(i as isize)).blk;
                *fresh83 =
                    rans_uncompress((*bc.offset(i as isize)).blk,
                                    (*bc.offset(i as isize)).sz,
                                    &mut (*bu.offset(i as isize)).sz as
                                        *mut uint32_t, order);
                i += 1
            }
            gettimeofday(&mut tv3 as *mut timeval, 0 as *mut timezone);
            i = 0i32;
            while i < nb {
                if (*b.offset(i as isize)).sz != (*bu.offset(i as isize)).sz
                       ||
                       0 !=
                           memcmp((*b.offset(i as isize)).blk as
                                      *const c_void,
                                  (*bu.offset(i as isize)).blk as
                                      *const c_void,
                                  (*b.offset(i as isize)).sz as c_ulong)
                {
                    panic!("Abort");
                }
                free((*bc.offset(i as isize)).blk as *mut c_void);
                free((*bu.offset(i as isize)).blk as *mut c_void);
                i += 1
            }
            free(bc as *mut c_void);
            free(bu as *mut c_void);
            fprintf(stderr,
                    (*::std::mem::transmute::<&[u8; 58],
                                              &mut [c_char; 58]>(b"%5.1f MB/s enc, %5.1f MB/s dec\t %lld bytes -> %lld bytes\n\x00")).as_mut_ptr(),
                    in_sz as c_double /
                        ((tv2.tv_sec - tv1.tv_sec) *
                             1000000i32 as c_long + tv2.tv_usec -
                             tv1.tv_usec) as c_double,
                    in_sz as c_double /
                        ((tv3.tv_sec - tv2.tv_sec) *
                             1000000i32 as c_long + tv3.tv_usec -
                             tv2.tv_usec) as c_double,
                    in_sz as c_longlong, out_sz as c_longlong);
             
        }
        exit(0i32);
    } else { return 0i32 };
}