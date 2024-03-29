#[cfg(feature = "log_color")]
#[macro_export]
macro_rules! log_info {
    ($msg:tt) => {
        println!(
            "\x1b[94m[{}]\x1b[m \x1b[36m[Thread ID: {}]\x1b[m \x1b[1m\x1b[32mNSP Info:\x1b[m {}",
            chrono::Local::now().format("%Y-%m-%d %H:%M:%S"),
            std::thread::current().id().as_u64(),
            $msg
        )
    };

    ($fmt_string:tt, $($args:tt)*) => {
        println!(
            concat!(
                "\x1b[94m[{}]\x1b[m \x1b[36m[Thread ID: {}]\x1b[m \x1b[1m\x1b[32mNSP Info:\x1b[m ",
                $fmt_string
            ),
            chrono::Local::now().format("%Y-%m-%d %H:%M:%S"),
            std::thread::current().id().as_u64(),
            $($args)*
        )
    };
}

pub fn log_info(msg: &str) {
    log_info!(msg);
}

#[cfg(feature = "log_color")]
#[macro_export]
macro_rules! log_warning {
    ($msg:tt) => {
        println!(
            "\x1b[94m[{}]\x1b[m \x1b[36m[Thread ID: {}]\x1b[m \x1b[1m\x1b[33mNSP Warning:\x1b[m {}",
            chrono::Local::now().format("%Y-%m-%d %H:%M:%S"),
            std::thread::current().id().as_u64(),
            $msg
        )
    };

    ($fmt_string:tt, $($args:tt)*) => {
        println!(
            concat!(
                "\x1b[94m[{}]\x1b[m \x1b[36m[Thread ID: {}]\x1b[m \x1b[1m\x1b[33mNSP Warning:\x1b[m ",
                $fmt_string
            ),
            chrono::Local::now().format("%Y-%m-%d %H:%M:%S"),
            std::thread::current().id().as_u64(),
            $($args)*
        )
    };
}

#[cfg(feature = "log_color")]
#[macro_export]
macro_rules! log_error {
    ($msg:tt) => {
        println!(
            "\x1b[94m[{}]\x1b[m \x1b[36m[Thread ID: {}]\x1b[m \x1b[1m\x1b[31mNSP Error:\x1b[m {}",
            chrono::Local::now().format("%Y-%m-%d %H:%M:%S"),
            std::thread::current().id().as_u64(),
            $msg
        )
    };

    ($fmt_string:tt, $($args:tt)*) => {
        println!(
            concat!(
                "\x1b[94m[{}]\x1b[m \x1b[36m[Thread ID: {}]\x1b[m \x1b[1m\x1b[31mNSP Error:\x1b[m ",
                $fmt_string
            ),
            chrono::Local::now().format("%Y-%m-%d %H:%M:%S"),
            std::thread::current().id().as_u64(),
            $($args)*
        )
    };
}

#[cfg(not(feature = "log_color"))]
#[macro_export]
macro_rules! log_info {
    ($msg:tt) => {
        println!(
            "[{}] [Thread ID: {}] NSP Info: {}",
            chrono::Local::now().format("%Y-%m-%d %H:%M:%S"),
            std::thread::current().id().as_u64(),
            $msg
        )
    };

    ($fmt_string:tt, $($args:tt)*) => {
        println!(
            concat!(
                "[{}] [Thread ID: {}] NSP Info: ",
                $fmt_string
            ),
            chrono::Local::now().format("%Y-%m-%d %H:%M:%S"),
            std::thread::current().id().as_u64(),
            $($args)*
        )
    };
}

#[cfg(not(feature = "log_color"))]
#[macro_export]
macro_rules! log_warning {
    ($msg:tt) => {
        println!(
            "[{}] [Thread ID: {}] NSP Warning: {}",
            chrono::Local::now().format("%Y-%m-%d %H:%M:%S"),
            std::thread::current().id().as_u64(),
            $msg
        )
    };

    ($fmt_string:tt, $($args:tt)*) => {
        println!(
            concat!(
                "[{}] [Thread ID: {}] NSP Warning: ",
                $fmt_string
            ),
            chrono::Local::now().format("%Y-%m-%d %H:%M:%S"),
            std::thread::current().id().as_u64(),
            $($args)*
        )
    };
}

#[cfg(not(feature = "log_color"))]
#[macro_export]
macro_rules! log_error {
    ($msg:tt) => {
        println!(
            "[{}] \x1b[36m[Thread ID: {}] \x1b[31mNSP Error: {}",
            chrono::Local::now().format("%Y-%m-%d %H:%M:%S"),
            std::thread::current().id().as_u64(),
            $msg
        )
    };

    ($fmt_string:tt, $($args:tt)*) => {
        println!(
            concat!(
                "[{}] \x1b[36m[Thread ID: {}] \x1b[31mNSP Error: ",
                $fmt_string
            ),
            chrono::Local::now().format("%Y-%m-%d %H:%M:%S"),
            std::thread::current().id().as_u64(),
            $($args)*
        )
    };
}

#[cfg(feature = "log_color")]
#[macro_export]
macro_rules! value_or_error {
    ($exp:expr, $err_value:expr, $msg:tt) => {
        match $exp {
            Err(e) => {
                println!(
                    "\x1b[94m[{}]\x1b[m \x1b[36m[Thread ID: {}]\x1b[m \x1b[1m\x1b[31mNSP Error:\x1b[m {}: {}",
                    chrono::Local::now().format("%Y-%m-%d %H:%M:%S"),
                    std::thread::current().id().as_u64(),
                    $msg,
                    e
                );
                return Err($err_value);
            },
            Ok(o) => o,
        }
    };
}

#[cfg(feature = "log_color")]
#[macro_export]
macro_rules! opt_or_error {
    ($exp:expr, $err_value:expr, $msg:tt) => {
        match $exp {
            None => {
                println!(
                    "\x1b[94m[{}]\x1b[m \x1b[36m[Thread ID: {}]\x1b[m \x1b[1m\x1b[31mNSP Error:\x1b[m {}",
                    chrono::Local::now().format("%Y-%m-%d %H:%M:%S"),
                    std::thread::current().id().as_u64(),
                    $msg
                );
                return Err($err_value);
            },
            Some(o) => o,
        }
    };
}

#[cfg(not(feature = "log_color"))]
#[macro_export]
macro_rules! value_or_error {
    ($exp:expr, $err_value:expr, $msg:tt) => {
        match $exp {
            Err(e) => {
                println!(
                    "[{}] [Thread ID: {}] NSP Error: {}: {}",
                    chrono::Local::now().format("%Y-%m-%d %H:%M:%S"),
                    std::thread::current().id().as_u64(),
                    $msg,
                    e
                );
                return Err($err_value);
            }
            Ok(o) => o,
        }
    };
}
