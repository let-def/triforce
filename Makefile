# SPDX-License-Identifier: GPL-2.0-or-later
# Copyright (C) 2024 James Calligeros <jcalligeros99@gmail.com>

PREFIX ?= /usr/local
LIBDIR ?= $(PREFIX)/lib64

test-release:
	cargo build --release
	rm -f output.wav
	LV2_PATH=$(PWD)/lv2-release ffmpeg -i ~/norm-f32.wav -af 'lv2=p=https\\://chadmed.au/triforce'  output.wav

test-debug:
	cargo build
	rm -f output.wav
	LV2_PATH=$(PWD)/lv2-debug ffmpeg -i ~/norm-f32.wav -af 'lv2=p=https\\://chadmed.au/triforce'  output.wav

default:
	cargo build --release

install:
	install -dDm0755 $(DESTDIR)/$(LIBDIR)/lv2/triforce.lv2/
	install -pm0755 target/release/libtriforce.so $(DESTDIR)/$(LIBDIR)/lv2/triforce.lv2/triforce.so
	install -pm0644 triforce.ttl $(DESTDIR)/$(LIBDIR)/lv2/triforce.lv2/triforce.ttl
	install -pm0644 manifest.ttl $(DESTDIR)/$(LIBDIR)/lv2/triforce.lv2/manifest.ttl

uninstall:
	rm -rf $(DESTDIR)/$(LIBDIR)/lv2/triforce.lv2/
