#![allow ( dead_code )]
#![allow ( mutable_transmutes )]
#![allow ( non_camel_case_types )]
#![allow ( non_snake_case )]
#![allow ( non_upper_case_globals )]
#![allow ( unused_mut )]
#![feature ( extern_types )]
#![feature ( i128_type )]
#![feature ( libc )]
#![feature ( offset_to )]
#![feature(start)]
extern crate brans;
use brans::rans_static::rans_compress;
use brans::rans_static::rans_uncompress;
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





pub fn main() {
    unsafe{main_0()};
}
pub unsafe extern "C" fn main_0() -> c_int {
    let mut opt: c_int = 0;
    let mut order: c_int = 1i32;
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
        len = 10000000i32 as size_t;
        let ref mut fresh80 = (*b.offset(nb as isize)).blk;
        *fresh80 = malloc(len) as *mut c_uchar;
        (*b.offset(nb as isize)).sz = len as uint32_t;
        xlen = len;
        let mut rnd_state = 1337u32;
        while 0 != xlen {
            rnd_state = rnd_state.wrapping_add(120720357);
            rnd_state %= 992687;
            xlen = xlen.wrapping_sub(1);
            *(*b.offset(nb as isize)).blk.offset(xlen as isize) =
                (rnd_state  as c_ulong & 63i32 as c_ulong) as c_uchar
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
