$(function () {
  var sidebarWidth = 260;
  
  function updateReportBodyPadding() {
    var toggleNavWidth = $(".toggleNav").outerWidth();
    var totalPadding = sidebarWidth + toggleNavWidth;
    $("#report_body").css("padding-left", totalPadding * 1.02 + "px");
    $(".toggleNav").css("left", sidebarWidth + "px");
  }
  
  $(window).load(function() {
    updateReportBodyPadding();
  });
  
  $(".fold1").click(function () {
    $(".jumpto-subnav").animate({ marginLeft: "-" + sidebarWidth + "px" });
    $(".toggleNav").animate({ left: "20px" });
    $(this).addClass("close");
    $(".fold2").removeClass("close");
  });
  $(".fold2").click(function () {
    $(".jumpto-subnav").animate({ marginLeft: "0px" });
    $(".toggleNav").animate({ left: sidebarWidth + "px" });
    $(this).addClass("close");
    $(".fold1").removeClass("close");
  });
});
