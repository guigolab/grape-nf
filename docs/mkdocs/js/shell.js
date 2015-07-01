$(".shell").each(function(index) {
  lines = $(this).text().trim().split("\n")
  $.each(lines, function(line) {
      lines[line]="<span class='line'>" + lines[line] + "</span>";
  });
  
  $(this).html(lines.join("\n"));
});