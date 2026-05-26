# Tier policy

## Tiers

msolve provides three tiers of target support:

- msolve provides no guarantees about tier 3 targets; they exist in the codebase, but may or may not build.
- msolve's continuous integration checks that tier 2 targets will always build, but they may or may not pass tests.
- msolve's continuous integration checks that tier 1 targets will always build and pass tests.

Adding a new tier 3 target imposes minimal requirements; we focus primarily on avoiding disruption to other ongoing msolve development.

Tier 2 and tier 1 targets place work on msolve project developers as a whole, to avoid breaking the target. Thus, these tiers require commensurate and ongoing efforts from the maintainers of the target, to demonstrate value and to minimize any disruptions to ongoing msolve development.

This policy defines the requirements for accepting a proposed target at a given level of support.

## Software platforms

### Tier 1
- x86_64 Linux
- x86_64 macOS
- ARM64 macOS

### Tier 2
- ARM64 Linux
- ARMv6 Linux
- ARMv7 Linux
- i686 Linux
- PPC64le Linux
- RISCV64 Linux
- x86_64 FreeBSD
- ARM64 FreeBSD

### Tier 3
- ARM64 Windows
- x86_64 Windows
- ARM64 Android
