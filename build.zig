const Builder = @import("std").build.Builder;

pub fn build(b: *Builder) void {
    const lib_cflags = [_][]const u8{"-std=c99"};

    const exe = b.addExecutable("zig-tracer", "src/main.zig");
    exe.setBuildMode(b.standardReleaseOptions());
    exe.addCSourceFile("src/pixel.c", lib_cflags);
    exe.addIncludeDir("src/");
    exe.linkSystemLibrary("c");
    b.default_step.dependOn(&exe.step);
    b.installArtifact(exe);
    if (exe.target.isDarwin()) {
        exe.addIncludeDir("/usr/local/Cellar/sdl2/2.0.9_1/include/SDL2");
    } else {
        //assuming this means linux
        exe.addIncludeDir("/usr/include/SDL2");
    }
    exe.linkSystemLibrary("SDL2");

    const run = b.step("run", "Run the project");
    const run_cmd = exe.run();
    run.dependOn(&run_cmd.step);
}
