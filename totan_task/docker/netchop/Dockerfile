# NetChop 3.1 (IEDB_NetChop-3.0 distribution) in a reproducible Linux container.
#
# Why a container: the shipped `netChop` is an ELF 32-bit i386 statically linked
# Linux executable. It cannot run natively on Windows or on Apple silicon, so a
# pinned Linux image is the only way to make the Week 2 processing run
# reproducible across machines.
#
# LICENSING: the NetChop distribution is NOT redistributed in this repository.
# It is licensed software (NetMHC 3.0 terms + Non-Profit OSL 3.0) and this repo
# is public. The tarball must be obtained separately and placed in the build
# context. See docs/week2_processing_reference_decision.md.
#
# Build (from repo root, with the tarball copied into docker/netchop/):
#   docker build -t sparq/netchop:3.1 docker/netchop
#
# Verify the install (the vendor's own fixture is unusable - see entrypoint.sh):
#   docker run --rm sparq/netchop:3.1 selftest
#
# Predict (20S network, -v 1):
#   docker run --rm -v "$PWD:/work" sparq/netchop:3.1 \
#       netChop -t 0.500000 -v 1 /work/input.fsa

# shell: process substitution needs bash, selftest uses it
FROM debian:bookworm-20241111-slim

# tcsh and gawk are required by the vendor scripts; python3 drives our wrapper.
# The netChop binary itself is statically linked, so no 32-bit userland libs are
# needed - only kernel ia32 support, which x86_64 Linux provides by default.
RUN apt-get update \
 && apt-get install -y --no-install-recommends \
        tcsh \
        bash \
        gawk \
        python3 \
        python3-pip \
        ca-certificates \
 && rm -rf /var/lib/apt/lists/*

ARG NETCHOP_TARBALL=IEDB_NetChop-3.0.tar.gz
WORKDIR /opt

COPY ${NETCHOP_TARBALL} /tmp/netchop.tar.gz
RUN tar -xzf /tmp/netchop.tar.gz -C /opt \
 && rm /tmp/netchop.tar.gz

# The vendor layout is /opt/netchop/netchop_3_1_executable/{bin,etc,data,test}.
# NETCHOP must point at that directory; TMPDIR needs the sticky bit.
ENV NETCHOP=/opt/netchop/netchop_3_1_executable \
    TMPDIR=/tmp \
    PATH=/opt/netchop/netchop_3_1_executable/bin:$PATH

RUN chmod 1777 /tmp \
 && chmod -R a+rx /opt/netchop/netchop_3_1_executable/bin

COPY entrypoint.sh /usr/local/bin/entrypoint.sh
RUN chmod +x /usr/local/bin/entrypoint.sh

WORKDIR /work
ENTRYPOINT ["/usr/local/bin/entrypoint.sh"]
CMD ["selftest"]
