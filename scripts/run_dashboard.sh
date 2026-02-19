#!/bin/sh
docker run --rm --name aspire \
  -p 18888:18888 \
  -p 4317:18889 \
  mcr.microsoft.com/dotnet/aspire-dashboard:latest
