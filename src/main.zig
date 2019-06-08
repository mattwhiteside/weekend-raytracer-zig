const Camera = @import("camera.zig").Camera;
const hitable = @import("hitable.zig");
const mat = @import("material.zig");
const Material = mat.Material;
const Lambertian = mat.Lambertian;
const Metal = mat.Metal;
const std = @import("std");
const debug = std.debug;
const ArrayList = std.ArrayList;
const rand = std.rand;
const Ray = @import("ray.zig").Ray;
const Vec3f = @import("vector.zig").Vec3f;

const HitRecord = hitable.HitRecord;
const Sphere = hitable.Sphere;
const World = hitable.World;

const c = @cImport({
    @cInclude("SDL.h");
});

// See https://github.com/zig-lang/zig/issues/565
// SDL_video.h:#define SDL_WINDOWPOS_UNDEFINED         SDL_WINDOWPOS_UNDEFINED_DISPLAY(0)
// SDL_video.h:#define SDL_WINDOWPOS_UNDEFINED_DISPLAY(X)  (SDL_WINDOWPOS_UNDEFINED_MASK|(X))
// SDL_video.h:#define SDL_WINDOWPOS_UNDEFINED_MASK    0x1FFF0000u
const SDL_WINDOWPOS_UNDEFINED = @bitCast(c_int, c.SDL_WINDOWPOS_UNDEFINED_MASK);

const window_width: c_int = 640;
const window_height: c_int = 320;
const num_samples: i32 = 256;
const max_depth: i32 = 32;

// For some reason, this isn't parsed automatically. According to SDL docs, the
// surface pointer returned is optional!
extern fn SDL_GetWindowSurface(window: *c.SDL_Window) ?*c.SDL_Surface;
extern fn setPixel(surf: *c.SDL_Surface, x: c_int, y: c_int, pixel: u32) void;

fn colorNormal(r: Ray, w: *const World) Vec3f {
    const maybe_hit = w.hit(r, 0.001, 10000.0);
    if (maybe_hit) |hit| {
        const n = hit.n.makeUnitVector();
        return n.add(Vec3f.one()).mul(0.5);
    } else {
        const unit_direction = r.direction.makeUnitVector();
        const t = 0.5 * (unit_direction.y + 1.0);
        return Vec3f.new(1.0, 1.0, 1.0).mul(1.0 - t).add(Vec3f.new(0.5, 0.7, 1.0).mul(t));
    }
}

fn colorAlbedo(r: Ray, w: *const World) Vec3f {
    const maybe_hit = w.hit(r, 0.001, 10000.0);
    if (maybe_hit) |hit| {
        return switch (hit.material) {
            Material.Lambertian => |l| l.albedo,
            Material.Metal => |m| m.albedo,
            Material.Dielectric => |l| Vec3f.one(),
        };
    } else {
        const unit_direction = r.direction.makeUnitVector();
        const t = 0.5 * (unit_direction.y + 1.0);
        return Vec3f.new(1.0, 1.0, 1.0).mul(1.0 - t).add(Vec3f.new(0.5, 0.7, 1.0).mul(t));
    }
}

fn colorDepthHelper(r: Ray, w: *const World, random: *rand.Random, depth: i32) i32 {
    const maybe_hit = w.hit(r, 0.001, 10000.0);
    if (maybe_hit) |hit| {
        if (depth < max_depth) {
            const scatter = switch (hit.material) {
                Material.Lambertian => |l| l.scatter(hit, random),
                Material.Metal => |m| m.scatter(r, hit, random),
                Material.Dielectric => |d| d.scatter(r, hit, random),
            };
            return colorDepthHelper(scatter.ray, w, random, depth + 1);
        } else {
            return depth; // reached max depth
        }
    } else {
        return depth; // hit the sky
    }
}

fn colorDepth(r: Ray, w: *const World, random: *rand.Random) Vec3f {
    const depth = colorDepthHelper(r, w, random, 0);
    return Vec3f.new(@intToFloat(f32, depth) / @intToFloat(f32, max_depth), 0.0, 0.0);
}

fn colorScattering(r: Ray, w: *const World, random: *rand.Random) Vec3f {
    const maybe_hit = w.hit(r, 0.001, 10000.0);
    if (maybe_hit) |hit| {
        const scatter = switch (hit.material) {
            Material.Lambertian => |l| l.scatter(hit, random),
            Material.Metal => |m| m.scatter(r, hit, random),
            Material.Dielectric => |d| d.scatter(r, hit, random),
        };
        const dir = scatter.ray.direction.makeUnitVector();
        return dir.add(Vec3f.one()).mul(0.5);
    } else {
        const unit_direction = r.direction.makeUnitVector();
        const t = 0.5 * (unit_direction.y + 1.0);
        return Vec3f.new(1.0, 1.0, 1.0).mul(1.0 - t).add(Vec3f.new(0.5, 0.7, 1.0).mul(t));
    }
}

fn color(r: Ray, world: *const World, random: *rand.Random, depth: i32) Vec3f {
    const maybe_hit = world.hit(r, 0.001, 10000.0);
    if (maybe_hit) |hit| {
        if (depth < max_depth) {
            const scatter = switch (hit.material) {
                Material.Lambertian => |l| l.scatter(hit, random),
                Material.Metal => |m| m.scatter(r, hit, random),
                Material.Dielectric => |d| d.scatter(r, hit, random),
            };
            return color(scatter.ray, world, random, depth + 1).elementwiseMul(scatter.attenuation);
        } else {
            return Vec3f.zero();
        }
    } else {
        const unit_direction = r.direction.makeUnitVector();
        const t = 0.5 * (unit_direction.y + 1.0);
        return Vec3f.new(1.0, 1.0, 1.0).mul(1.0 - t).add(Vec3f.new(0.5, 0.7, 1.0).mul(t));
    }
}

fn toBgra(r: u32, g: u32, b: u32) u32 {
    return 255 << 24 | r << 16 | g << 8 | b;
}

pub fn main() !void {
    if (c.SDL_Init(c.SDL_INIT_VIDEO) != 0) {
        c.SDL_Log(c"Unable to initialize SDL: %s", c.SDL_GetError());
        return error.SDLInitializationFailed;
    }
    defer c.SDL_Quit();

    const window = c.SDL_CreateWindow(c"weekend raytracer", SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED, window_width, window_height, c.SDL_WINDOW_OPENGL) orelse {
        c.SDL_Log(c"Unable to create window: %s", c.SDL_GetError());
        return error.SDLInitializationFailed;
    };
    defer c.SDL_DestroyWindow(window);

    const surface = SDL_GetWindowSurface(window) orelse {
        c.SDL_Log(c"Unable to get window surface: %s", c.SDL_GetError());
        return error.SDLInitializationFailed;
    };

    // Ray tracing takes place here

    const lookfrom = Vec3f.new(16.0, 2.0, 4.0);
    const lookat = Vec3f.new(0.0, 0.0, 0.0);
    const vfov = 15.0;
    const focus_distance = lookfrom.sub(lookat).length();
    const aperture = 0.2;
    // 640 by 320
    const aspect_ratio = @intToFloat(f32, window_width) / @intToFloat(f32, window_height);
    const camera = Camera.new(lookfrom, lookat, Vec3f.new(0.0, 1.0, 0.0), vfov, aspect_ratio, aperture, focus_distance);

    var world = World.init();
    defer world.deinit();

    try world.spheres.append(Sphere.new(Vec3f.new(0.0, -1000.0, -1.0), 1000.0, Material.lambertian(Vec3f.new(0.5, 0.5, 0.5))));
    try world.spheres.append(Sphere.new(Vec3f.new(0.0, 1.0, 0.0), 1.0, Material.dielectric(1.5)));
    try world.spheres.append(Sphere.new(Vec3f.new(-4.0, 1.0, 0.0), 1.0, Material.lambertian(Vec3f.new(0.4, 0.2, 0.1))));
    try world.spheres.append(Sphere.new(Vec3f.new(4.0, 1.0, 0.0), 1.0, Material.metal(Vec3f.new(0.7, 0.6, 0.5), 0.0)));

    var prng = rand.DefaultPrng.init(0);

    const sphere_offset = Vec3f.new(4.0, 0.2, 0.0);
    var i: i32 = -5;
    while (i < 5) : (i += 1) {
        var j: i32 = -5;
        while (j < 5) : (j += 1) {
            const a = @intToFloat(f32, i);
            const b = @intToFloat(f32, j);
            const center = Vec3f.new(a + 0.9 * prng.random.float(f32), 0.2, b + 0.9 * prng.random.float(f32));
            const choose_mat = prng.random.float(f32);
            if (center.sub(sphere_offset).length() > 0.9) {
                if (choose_mat < 0.8) {
                    // diffuse
                    const random_albedo = Vec3f.new(prng.random.float(f32), prng.random.float(f32), prng.random.float(f32));
                    try world.spheres.append(Sphere.new(center, 0.2, Material.lambertian(random_albedo)));
                } else if (choose_mat < 0.95) {
                    // metal
                    const random_albedo = Vec3f.new(prng.random.float(f32), prng.random.float(f32), prng.random.float(f32));
                    try world.spheres.append(Sphere.new(center, 0.2, Material.metal(random_albedo, 0.5 * prng.random.float(f32))));
                } else {
                    try world.spheres.append(Sphere.new(center, 0.2, Material.dielectric(1.5)));
                }
            }
        }
    }

    {
        _ = c.SDL_LockSurface(surface);
        var idx: i32 = 0;
        while (idx < window_width * window_height) : (idx += 1) {
            const w = @mod(idx, window_width);
            const h = @divTrunc(idx, window_width);
            var sample: i32 = 0;
            var color_accum = Vec3f.zero();

            while (sample < num_samples) : (sample += 1) {
                const u = (@intToFloat(f32, w) + prng.random.float(f32)) / @intToFloat(f32, window_width);
                const v = (@intToFloat(f32, h) + prng.random.float(f32)) / @intToFloat(f32, window_height);

                const r = camera.makeRay(&prng.random, u, v);
                const color_sample = color(r, &world, &prng.random, 0);
                // const color_sample = colorScattering(r, &world, &prng.random);
                // const color_sample = colorDepth(r, &world, &prng.random);
                // const color_sample = colorNormal(r, &world);
                // const color_sample = colorAlbedo(r, &world);
                color_accum = color_accum.add(color_sample);
            }
            color_accum = color_accum.mul(1.0 / @intToFloat(f32, num_samples));
            setPixel(surface, w, window_height - h - 1, toBgra(@floatToInt(u32, 255.99 * color_accum.x), @floatToInt(u32, 255.99 * color_accum.y), @floatToInt(u32, 255.99 * color_accum.z)));
        }
        c.SDL_UnlockSurface(surface);
    }

    if (c.SDL_UpdateWindowSurface(window) != 0) {
        c.SDL_Log(c"Error updating window surface: %s", c.SDL_GetError());
        return error.SDLUpdateWindowFailed;
    }

    var running = true;
    while (running) {
        var event: c.SDL_Event = undefined;
        while (c.SDL_PollEvent(&event) != 0) {
            switch (event.@"type") {
                c.SDL_QUIT => {
                    running = false;
                },
                else => {},
            }
        }

        c.SDL_Delay(16);
    }
}
