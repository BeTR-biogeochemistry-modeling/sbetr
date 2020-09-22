#ifdef NDEBUG
#define SHR_ASSERT(assert, msg, bstatus)
#define SHR_ASSERT_ALL(assert, msg, bstatus)
#define SHR_ASSERT_ANY(assert, msg, bstatus)
#define SHR_ASSERT_ALL_EXT(assert, msg)
#else
#define SHR_ASSERT(assert, msg, bstatus) call shr_assert(assert, msg, bstatus); if(bstatus%check_status())return
#define SHR_ASSERT_ALL(assert, msg, bstatus) call shr_assert_all(assert, msg, bstatus); if(bstatus%check_status())return
#define SHR_ASSERT_ALL_EXT(assert, msg) call shr_assert_all_ext(assert, msg)
#define SHR_ASSERT_ANY(assert, msg, bstatus) call shr_assert_any(assert, msg, bstatus); if(bstatus%check_status())return
#endif
