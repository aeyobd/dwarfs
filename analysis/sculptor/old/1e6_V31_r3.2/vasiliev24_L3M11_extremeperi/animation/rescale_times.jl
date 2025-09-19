using LilGuys

out = Output("..")
times = out.times
times_new = times * T2GYR / 0.97779

open("times.txt", "w") do f
    for time in times_new
        println(f, time)
    end
end
