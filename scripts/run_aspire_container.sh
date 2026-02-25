#!/bin/sh
docker run --rm --name aspire \
  -p 18888:18888 \
  -p 18891:18891 \
  -p 18889:18889 \
  -e ASPIRE_DASHBOARD_UNSECURED_ALLOW_ANONYMOUS=true \
  -e ASPIRE_ALLOW_UNSECURED_TRANSPORT=true \
  -e ASPIRE_DASHBOARD_MCP_ENDPOINT_URL="http://0.0.0.0:18891" \
  -e Dashboard__Mcp__AuthMode=Unsecured \
  -e Dashboard__Frontend__EndpointUrls=http://localhost:18888 \
  -e Dashboard__Frontend__PublicUrl=http://localhost:18888 \
  mcr.microsoft.com/dotnet/aspire-dashboard:latest
