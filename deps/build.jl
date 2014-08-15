using BinDeps

@BinDeps.setup

libwsmp = library_dependency("libwsmp")

libdir = joinpath(BinDeps.depsdir(libwsmp),"usr","lib")

provides(SimpleBuild,
         (@build_steps begin
             ChangeDirectory(libdir)
             `make`
         end),libwsmp,os = :Unix)
         
@BinDeps.install [:libwsmp => :libwsmp]
