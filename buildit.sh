#!/usr/bin/env bash
eupspkg -e -v 1 fetch
eupspkg -e -v 1 prep
eupspkg -e -v 1 config
eupspkg -e -v 1 build
eupspkg -e -v 1 install
