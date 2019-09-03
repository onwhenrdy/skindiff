library(dskin)

# setup
base.dir <- getwd()
project.name <- "test"
working.dir <- project.name
if (!dir.exists(working.dir))
{
  dir.create(working.dir, recursive = TRUE)
}

# CONFIG parameters
p <- dskin.parameter(project.name,
                    n.layers = 2,
                    layer.names = c("Stratum corneum", "Deeper skin layers"),
                    vehicle.name = "Vehicle",
                    sink.name = "Blood level")
p$sys$max_module = 200
p$sys$sim_time = 24*60*18
p$sys$disc_scheme = "BK"
p$sys$resolution = 4
p$log$scaling = "ng"
# vehicle
p$compartments$vehicle$c_init = 127.2727
p$compartments$vehicle$h = 110
p$compartments$vehicle$D = 9.266667
p$compartments$vehicle$app_area = 15.0
# SC
p$compartments$layers$h[1] = 190
p$compartments$layers$D[1] = 28.2539
p$compartments$layers$K[1] = 421.543
p$compartments$layers$cross_section[1] = 0.001
# DSL
p$compartments$layers$h[2] = 200
p$compartments$layers$D[2] = 5767.783
p$compartments$layers$K[2] = 0.04719648
p$compartments$layers$cross_section[2] = 0.3
p$PK$enabled = TRUE
p$PK$t_half = 24
p$compartments$sink$Vd = 25 * 1000 * 75    # 20 - 25 l /kg ~ 25 * 1000 * 75
p$compartments$vehicle$remove_after = 24 * 60 * 16
p$compartments$vehicle$replace_after = 24 * 60

# RUNTIME CODE
setwd(working.dir)

# simulate
sim.res <- dskin.simulate(p)
cat("Runtime:", dskin.runtime.str(sim.res))

# plots
mass.file <- paste(project.name, "_masses.pdf", sep="")
pdf(mass.file)
dskin.plot.masses(p, sim.res, max.cols = 2, type="l", lwd=2)
dev.off()

conc.file <- paste(project.name, "_concs.pdf", sep="")
pdf(conc.file)
dskin.plot.concs(p, sim.res, max.cols = 2, type="l", lwd=2, max.profiles = 10)
dev.off()

# movies
m.file <- dskin.movie.masses(p, sim.res, max.cols = 2, type="l", lwd=7, col="orange", cex.axis = 2, cex.lab = 1.7, cex.main = 2)
cat("Movie produced:", m.file)
m.file.gif <- dskin.movie.to.gif(p, width = 800, fps = 10)
cat("Gif produced:", m.file.gif)

c.file <- dskin.movie.concs(p, sim.res, max.cols = 2, type="l", lwd = 2, col = "red")
cat("Movie produced:", c.file)
c.file.gif <- dskin.movie.to.gif(p, dskin.class = "conc", width = 720, fps = 10)
cat("Gif produced:", c.file.gif)

# return to base dir
setwd(base.dir)

